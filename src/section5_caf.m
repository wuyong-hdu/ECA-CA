function metric = section5_caf(ref, residual, outputStart, interval, cfg)
% CAF neighborhood centered at the configured physical target Doppler.
% Adjacent Doppler samples are spaced by 1/(actual CAF duration).
% Direct DFT evaluation keeps the center fixed as the CPI length changes.
% Delayed reference samples before the CPI are zero, matching the legacy CAF.
first = interval(1);
count = interval(2)-first+1;
L = cfg.cafAverageLength;
assert(mod(count,L)==0 && count>0);
assert(first>=1 && first==floor(first));
offset = first-outputStart+1;
assert(offset>=1 && offset+count-1<=numel(residual));
assert(interval(2)<=numel(ref));
signal = residual(offset:offset+count-1);
blocks = count/L;
gate = cfg.cafDopplerHalfWidth;
validateattributes(gate,{'numeric'},{'scalar','integer','positive'});
assert(blocks>=2*gate+1,'CAF duration is too short for the Doppler neighborhood.');
validateattributes(cfg.targetRangeBin,{'numeric'},{'scalar','integer','positive','<=',cfg.rangeBins});
spacing = cfg.fs/count;
frequencies = cfg.targetDopplerHz+(-gate:gate)*spacing;
assert(all(frequencies>=-cfg.fs/(2*L) & frequencies<cfg.fs/(2*L)), ...
    'The Doppler neighborhood must lie within the slow-time Nyquist interval.');
% Center time of each nonoverlapping average, relative to the evaluated interval.
times = ((0:blocks-1)'*L+(L-1)/2)/cfg.fs;
weights = exp(-1i*2*pi*(times*frequencies));
localCAF = zeros(2*gate+1,cfg.rangeBins);
% Native FFT grid for the overview map; keep the target-centered auxiliary grid
% unchanged. Save one bin beyond a requested limit when needed to cover it.
collectMap = isfield(cfg,'collectCAFMap') && cfg.collectCAFMap;
map = struct([]);
if collectMap
    limits = cfg.cafMapDopplerLimitsHz;
    validateattributes(limits,{'numeric'},{'vector','numel',2,'real','finite'});
    assert(limits(1)<limits(2));
    assert(limits(1)>=-cfg.fs/(2*L) && limits(2)<cfg.fs/(2*L));
    bins = max(-floor(blocks/2),floor(limits(1)/spacing)): ...
        min(ceil(blocks/2)-1,ceil(limits(2)/spacing));
    fftRows = bins+floor(blocks/2)+1;
    mapMagnitude = zeros(numel(bins),cfg.rangeBins);
    allMagnitude = zeros(blocks,cfg.rangeBins);
end
for bin=1:cfg.rangeBins
    delay = bin-1;
    product = complex(zeros(count,1));
    validFirst = max(first,delay+1);
    if validFirst<=interval(2)
        localFirst = validFirst-first+1;
        product(localFirst:end) = conj(ref(validFirst-delay:interval(2)-delay)) ...
            .*signal(localFirst:end);
    end
    slow = mean(reshape(product,L,blocks),1).';
    localCAF(:,bin) = abs(weights.'*slow)/blocks;
    if collectMap
        spectrum = fftshift(fft(slow));
        mapMagnitude(:,bin) = abs(spectrum(fftRows))/blocks;
        allMagnitude(:,bin) = abs(spectrum)/blocks;
    end
end
if collectMap
    map = struct('Magnitude',mapMagnitude,'RangeBins',1:cfg.rangeBins, ...
        'DopplerHz',bins*spacing,'DopplerBinSpacingHz',spacing, ...
        'RangeLimits',cfg.cafMapRangeLimits,'DopplerLimitsHz',limits, ...
        'Transform',"native_fft",'AllDopplerMagnitude',allMagnitude, ...
        'AllDopplerHz',(-floor(blocks/2):ceil(blocks/2)-1)*spacing);
end
center = gate+1;
slice = localCAF(center,:);
metric = struct('Slice',slice, ...
    'StartSample',first,'EndSample',interval(2),'DurationSeconds',count/cfg.fs, ...
    'DopplerHz',cfg.targetDopplerHz,'DopplerBinSpacingHz',spacing, ...
    'DopplerHalfWidth',gate, ...
    'LocalCAF',localCAF,'LocalDopplerHz',frequencies, ...
    'ReferenceBoundary',"zero_before_cpi",'RangeDopplerMap',map);
metric=section5_apply_ptr(metric,cfg);
end
