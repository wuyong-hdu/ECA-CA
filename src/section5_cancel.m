function [err, info] = section5_cancel(method, ref, sig, K, B, cfg, collectMemory, boundary)
% Clutter suppression kernels. No parallel workers or I/O.
% Memory is the high-water mark of named, explicitly allocated array payloads.
% N-sample inputs, final residual, temporaries and library workspace are excluded.
% FS-only boundary buffers and assembled extended-input copies are included.
% boundary contains RefBefore, RefAfter, SigBefore, SigAfter from the recording.
if nargin < 7, collectMemory = false; end
if nargin < 8, boundary = struct(); end
validateattributes(ref, {'double'}, {'column','finite','nonempty'});
validateattributes(sig, {'double'}, {'column','finite','numel',numel(ref)});
N = numel(ref);
validateattributes(K, {'numeric'}, {'scalar','integer','positive','<=',N});
validateattributes(B, {'numeric'}, {'scalar','integer','positive'});
assert(mod(N,B)==0, 'N must be divisible by B.');
switch string(method)
    case "ECA"
        [err,info] = full_eca(ref,sig,K,collectMemory);
    case "ECA-B"
        assert(N/B >= K, 'ECA-B requires at least K samples in each batch.');
        [err,info] = serial_ecab(ref,sig,K,B,collectMemory,cfg);
    case "ECA-CA"
        [err,info] = accumulated_eca(ref,sig,K,B,collectMemory);
    case "ECA-FS"
        [err,info] = frequency_eca(ref,sig,B,cfg,collectMemory,boundary);
    otherwise
        error('Unknown method: %s',method);
end
assert(all(isfinite(err)), 'Cancellation returned nonfinite samples.');
info.Method = string(method);
info.CPILength = N;
info.InputLength = N+info.PrefixSamples+info.SuffixSamples;
info.OutputLength = numel(err);
end

function info = initial_info()
info = struct('OutputStart',1,'InputStart',1,'PrefixSamples',0,'SuffixSamples',0, ...
    'BoundaryInputBytes',0,'PeakWorkingBytes',NaN,'PeakStage',"", ...
    'MinimumRcond',NaN,'FFTLength',NaN,'Segments',NaN,'FrequencyBatches',NaN, ...
    'PeakArrays',struct([]),'DelayBoundary',"full_cpi_circular", ...
    'CoefficientSolver',"backslash");
end

function info = account(info, vars, names, stage)
% whos reports payload bytes. Names list makes the inclusion convention explicit.
selected = vars(ismember({vars.name},names));
bytes = sum([selected.bytes]);
if isnan(info.PeakWorkingBytes) || bytes > info.PeakWorkingBytes
    info.PeakWorkingBytes = bytes;
    info.PeakStage = string(stage);
    info.PeakArrays = struct('Name',{selected.name},'Size',{selected.size}, ...
        'Class',{selected.class},'Bytes',{selected.bytes});
end
end

function [err,info] = full_eca(ref,sig,K,collect)
N = numel(ref);
info = initial_info();
names = {'V','G','rhs','w','clutter'};
V = complex(zeros(N,K));
for k=1:K
    V(k:N,k) = ref(1:N-k+1);
    if k>1, V(1:k-1,k) = ref(N-k+2:N); end
end
G = V' * V;
rhs = V' * sig;
w = G \ rhs;                     % Shared normal-equation solver for ECA/CA/B.
if collect
    info = account(info,whos,names,'solve');
    info.MinimumRcond = rcond(G);
end
clear G rhs
clutter = V*w;
if collect, info = account(info,whos,names,'reconstruction'); end
err = sig-clutter;
end

function [err,info] = serial_ecab(ref,sig,K,B,collect,cfg)
N = numel(ref); D = N/B;
info = initial_info(); info.Segments = B;
names = {'V','indices','referenceBlock','blockSig','G','rhs','w','clutter','residualBlock'};
implementation="continuous_history";
if isfield(cfg,'ecabImplementation'),implementation=string(cfg.ecabImplementation);end
assert(any(implementation==["legacy_batch","continuous_history"]));
if implementation=="legacy_batch"
    info.DelayBoundary="batch_circular";
end
err = complex(zeros(N,1));
minimumRcond = Inf;
for b=1:B
    first = (b-1)*D+1; last = b*D;
    V = complex(zeros(D,K));
    if implementation=="legacy_batch"
        referenceBlock=ref(first:last);
        for k=1:K
            V(k:D,k)=referenceBlock(1:D-k+1);
            if k>1,V(1:k-1,k)=referenceBlock(D-k+2:D);end
        end
    else
        for k=1:K
            indices = mod((first:last)-(k-1)-1,N)+1;
            V(:,k) = ref(indices);
        end
    end
    if collect, info = account(info,whos,names,'batch formation'); end
    clear indices referenceBlock
    blockSig = sig(first:last);
    G = V'*V;
    rhs = V'*blockSig;
    w=G\rhs;
    if collect
        info = account(info,whos,names,'batch solve');
        minimumRcond = min(minimumRcond,rcond(G));
    end
    clear G rhs
    clutter = V*w;
    residualBlock = blockSig-clutter;
    if collect, info = account(info,whos,names,'batch reconstruction'); end
    err(first:last) = residualBlock;
    clear V blockSig w clutter residualBlock
end
if collect, info.MinimumRcond = minimumRcond; end
end

function [err,info] = accumulated_eca(ref,sig,K,B,collect)
N = numel(ref); D = N/B; H = K-1;
L = D+H;                        % Exact block FFT length from Section 3.
info = initial_info(); info.FFTLength = L; info.Segments = B;
names = {'extended','history','offsets','A','currentSpectrum','conjugateA', ...
    'product','autoSpectrum','surBlocks','surSpectrum','crossSpectrum', ...
    'autoCorr','crossCorr','G','rhs','w','coefficientSpectrum','blocks','clutter'};
extended = complex(zeros(L,B));
extended(K:L,:) = reshape(ref,D,B);
if H>0
    offsets = (0:B-1)*D;
    history = mod((-H:-1).' + offsets,N)+1;
    extended(1:H,:) = ref(history);
    if collect, info = account(info,whos,names,'extended reference'); end
    clear history offsets
end
A = fft(extended,L,1);
extended(1:H,:) = 0;
currentSpectrum = fft(extended,L,1);
if collect, info = account(info,whos,names,'reference transforms'); end
clear extended
conjugateA = conj(A);
product = conjugateA .* currentSpectrum;
if collect, info = account(info,whos,names,'autocorrelation products'); end
clear currentSpectrum
autoSpectrum = sum(product,2);
clear product
surBlocks = complex(zeros(L,B));
surBlocks(K:L,:) = reshape(sig,D,B);
surSpectrum = fft(surBlocks,L,1);
if collect, info = account(info,whos,names,'surveillance transform'); end
clear surBlocks
product = conjugateA .* surSpectrum;
if collect, info = account(info,whos,names,'cross-correlation products'); end
clear conjugateA surSpectrum
crossSpectrum = sum(product,2);
clear product
autoCorr = ifft(autoSpectrum,L);
crossCorr = ifft(crossSpectrum,L);
autoCorr(1) = real(autoCorr(1));  % Zero-lag energy is real; remove FFT roundoff.
if collect, info = account(info,whos,names,'accumulated correlations'); end
clear autoSpectrum crossSpectrum
G = toeplitz(autoCorr(1:K),conj(autoCorr(1:K)));
rhs = crossCorr(1:K);
if collect, info = account(info,whos,names,'normal equations'); end
clear autoCorr crossCorr
w = G\rhs;
if collect
    info = account(info,whos,names,'global solve');
    info.MinimumRcond = rcond(G);
end
clear G rhs
coefficientSpectrum = fft(w,L,1);
product = A .* coefficientSpectrum;
if collect, info = account(info,whos,names,'reconstruction product'); end
clear A coefficientSpectrum
blocks = ifft(product,L,1);
if collect, info = account(info,whos,names,'reconstruction transform'); end
clear product
clutter = reshape(blocks(K:L,:),N,1);
if collect, info = account(info,whos,names,'valid clutter samples'); end
clear blocks
err = sig-clutter;
end

function [err,info] = frequency_eca(ref,sig,B,cfg,collect,boundary)
N = numel(ref); layout = section5_fs_layout(N,B,cfg);
L = layout.SegmentLength; step = layout.Step; M = layout.Segments;
assert(all(isfield(boundary,{'RefBefore','RefAfter','SigBefore','SigAfter'})), ...
    'ECA-FS requires real boundary samples; use section5_fs_layout to obtain their lengths.');
refBefore = boundary.RefBefore; refAfter = boundary.RefAfter;
sigBefore = boundary.SigBefore; sigAfter = boundary.SigAfter;
validateattributes(refBefore,{'double'},{'column','finite','numel',layout.PrefixSamples});
validateattributes(sigBefore,{'double'},{'column','finite','numel',layout.PrefixSamples});
validateattributes(refAfter,{'double'},{'column','finite','numel',layout.SuffixSamples});
validateattributes(sigAfter,{'double'},{'column','finite','numel',layout.SuffixSamples});
Bf = ceil(M/cfg.fsBatchTarget);
Mb = ceil(M/Bf);
assert(cfg.fsExtendedBatch>=Mb);
extension = floor((cfg.fsExtendedBatch-Mb)/2);
info = initial_info();
info.InputStart = layout.InputStart;
info.DelayBoundary="frequency_domain";info.CoefficientSolver="scalar_ls";
info.PrefixSamples = layout.PrefixSamples; info.SuffixSamples = layout.SuffixSamples;
info.SegmentStep = step; info.RetainedStart = layout.RetainedStart;
boundaryVars = whos('refBefore','refAfter','sigBefore','sigAfter');
info.BoundaryInputBytes = sum([boundaryVars.bytes]);
info.FFTLength = L; info.Segments = M; info.FrequencyBatches = Bf;
names = {'Rseg','Sseg','Rf','Sf','Ef','Eseg','Rk','Sk','Ek','Rp','Sp', ...
    'batchResult','validIndices','outIndices','refExtended','sigExtended', ...
    'refBefore','refAfter','sigBefore','sigAfter'};
% Assemble inside the timed call. Boundary input buffers remain live throughout;
% their aliases above are counted once, without counting the boundary struct again.
refExtended = [refBefore;ref;refAfter];
sigExtended = [sigBefore;sig;sigAfter];
Rseg = complex(zeros(L,M)); Sseg = complex(zeros(L,M));
for m=1:M
    first = (m-1)*step+1;
    Rseg(:,m) = refExtended(first:first+L-1);
    Sseg(:,m) = sigExtended(first:first+L-1);
end
if collect, info = account(info,whos,names,'extended input and segmentation'); end
clear refExtended sigExtended
Rseg(L-cfg.fsTailZeros+1:L,:) = 0;
Rf = fft(Rseg,[],1);
if collect, info = account(info,whos,names,'reference transform'); end
clear Rseg
Sf = fft(Sseg,[],1);
if collect, info = account(info,whos,names,'surveillance transform'); end
clear Sseg
Ef = complex(zeros(L,M));
% Process every frequency bin.
for k=1:L
    Rk = Rf(k,:).'; Sk = Sf(k,:).'; Ek = complex(zeros(M,1));
    for b=1:Bf
        first = (b-1)*Mb+1; last = min(b*Mb,M);
        left = max(1,first-extension); right = min(M,last+extension);
        Rp = Rk(left:right); Sp = Sk(left:right);
        denominator = real(Rp'*Rp)+cfg.fsLambda;
        if denominator==0, alpha=0; else, alpha=(Rp'*Sp)/denominator; end
        batchResult = Sk(first:last)-Rk(first:last)*alpha;
        Ek(first:last) = batchResult;
        if collect, info = account(info,whos,names,'frequency-bin estimation'); end
        clear Rp Sp batchResult
    end
    Ef(k,:) = Ek.';
    clear Rk Sk Ek
end
clear Rf Sf
Eseg = ifft(Ef,[],1);
if collect, info = account(info,whos,names,'inverse transforms'); end
clear Ef
validIndices = cfg.fsOutputStart:cfg.fsOutputStart+step-1;
err = complex(zeros(M*step,1));
for m=1:M
    outIndices = (m-1)*step+(1:step);
    err(outIndices) = Eseg(validIndices,m);
end
if collect, info = account(info,whos,names,'output assembly'); end
end
