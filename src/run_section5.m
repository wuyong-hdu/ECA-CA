function results = run_section5(cfg)
% RUN_SECTION5 Run the common experiment cases and save reproducible results.
%   results = run_section5();            % Full 14 cases, 1 warm-up + 50 timings.
%   cfg = section5_config(); ...; results = run_section5(cfg);
if nargin==0, cfg=section5_config(); end
validateattributes(cfg.repetitions,{'numeric'},{'scalar','integer','positive'});
validateattributes(cfg.warmups,{'numeric'},{'scalar','integer','nonnegative'});
validateattributes(cfg.startSample,{'numeric'},{'scalar','integer','positive'});
assert(cfg.fsLambda>=0);
validateattributes(cfg.ptrThresholdScale,{'numeric'},{'scalar','real','finite','positive'});
validateattributes(cfg.cafDopplerHalfWidth,{'numeric'},{'scalar','integer','positive'});
validateattributes(cfg.targetRangeBin,{'numeric'},{'scalar','integer','positive','<=',cfg.rangeBins});
[~,validTargetSlice]=section5_mask_caf(zeros(1,cfg.rangeBins),1:cfg.rangeBins,cfg.targetDopplerHz,cfg);
assert(validTargetSlice(cfg.targetRangeBin),'section5:MaskedTarget','The configured target lies inside the CAF mask.');
cases = cfg.cases;
assert(height(cases)>0 && height(unique(cases(:,{'N','K','B'})))==height(cases));
assert(all(cases.N>0 & cases.K>0 & cases.B>0));
assert(all(mod(cases.N,cases.B)==0 & cases.N./cases.B>=cases.K));
assert(all(mod(cases.N,4*cases.B)==0));

% N denotes the complete evaluated output length. Keep each method's native
% acquisition interval and pair its residual with the corresponding reference.
prefix=zeros(height(cases),1); suffix=prefix;
sampleIntervals=struct([]);
for c=1:height(cases)
    layout=section5_fs_layout(cases.N(c),cases.B(c),cfg);
    prefix(c)=layout.PrefixSamples; suffix(c)=layout.SuffixSamples;
    acquired=section5_input_intervals(cases.N(c),cases.B(c),cfg);
    if c==1,sampleIntervals=acquired;else,sampleIntervals(c)=acquired;end
end
readStart=min([cfg.startSample sampleIntervals.FSInputStart]);
readEnd=max([cfg.startSample+max(cases.N)-1 sampleIntervals.FSInputEnd]);
assert(readStart>=1,'The recording has insufficient real ECA-FS prefix samples.');
refData=read_channel(cfg.dataFile,cfg.referenceVariable,readStart,readEnd-readStart+1);
sigData=read_channel(cfg.dataFile,cfg.surveillanceVariable,readStart,readEnd-readStart+1);
assert(all(isfinite(refData)) && all(isfinite(sigData)));
cpiOffset=cfg.startSample-readStart+1;

% Each CAF uses all N samples, independent of other cases in this run.
lengths = unique(cases.N);
assert(all(mod(lengths,cfg.cafAverageLength)==0), ...
    'A full CPI must contain an integer number of CAF averaging blocks.');
assert(all(lengths>=(2*cfg.cafDopplerHalfWidth+1)*cfg.cafAverageLength));
intervals = [lengths ones(numel(lengths),1) lengths];

runName = sprintf('%s_%s',cfg.label,char(datetime('now','Format','yyyyMMdd_HHmmss_SSS')));
outDir = fullfile(cfg.outputRoot,runName);
assert(~isfolder(outDir),'Result directory already exists.');
mkdir(outDir);
results = struct('Config',cfg,'Cases',cases,'Intervals',intervals, ...
    'Metrics',table(),'Slices',{{}},'Raw',struct([]),'RunInfo',struct([]), ...
    'Complete',false,'OutputDirectory',outDir);
results.Environment = struct('MATLAB',version,'Computer',computer, ...
    'LibraryThreads',maxNumCompThreads,'Created',char(datetime('now')));
results.DataInterval=struct('ReadStart',readStart,'ReadEnd',readEnd, ...
    'CPIStart',cfg.startSample,'Convention',sampleIntervals(1).Convention, ...
    'Cases',sampleIntervals);
% Keep the exact code with the results so later source edits do not change provenance.
codeDir = fileparts(mfilename('fullpath'));
mkdir(fullfile(outDir,'source'));
sourceFiles = dir(fullfile(codeDir,'*.m'));
for f=1:numel(sourceFiles)
    copyfile(fullfile(codeDir,sourceFiles(f).name),fullfile(outDir,'source'));
end

cache = containers.Map('KeyType','char','ValueType','any');
methods = ["ECA","ECA-B","ECA-FS","ECA-CA"];
rows = struct([]);
try
    for c=1:height(cases)
        N=cases.N(c); K=cases.K(c); B=cases.B(c);
        cpiLast=cpiOffset+N-1;
        ref=refData(cpiOffset:cpiLast); sig=sigData(cpiOffset:cpiLast);
        intRow=find(lengths==N,1); interval=intervals(intRow,2:3);
        fprintf('Case %d/%d: N=%d, K=%d, B=%d\n',c,height(cases),N,K,B);
        if numel(results.Raw)<intRow || isempty(results.Raw(intRow).N)
            raw=section5_caf(ref,sig,1,interval,cfg);
            results.Raw(intRow).N=N;
            results.Raw(intRow).Metric=raw;
        end
        for m=1:numel(methods)
            method=methods(m);
            methodFirst=cfg.startSample;
            inputConvention="base_cpi"; sourceFrameIndex=NaN;
            if method=="ECA-FS"
                methodFirst=sampleIntervals(c).FSOutputStart;
                inputConvention=sampleIntervals(c).Convention;
                sourceFrameIndex=sampleIntervals(c).FrameIndex;
            end
            methodOffset=methodFirst-readStart+1; methodLast=methodOffset+N-1;
            methodRef=refData(methodOffset:methodLast);
            methodSig=sigData(methodOffset:methodLast);
            switch method
                case "ECA", key=sprintf('ECA_N%d_K%d',N,K);
                case "ECA-FS", key=sprintf('FS_N%d_B%d',N,B);
                otherwise, key=sprintf('%s_N%d_K%d_B%d',method,N,K,B);
            end
            reused=isKey(cache,key);
            if reused
                run=cache(key);
            else
                boundary=struct();
                if method=="ECA-FS"
                    boundary=struct( ...
                        'RefBefore',refData(methodOffset-prefix(c):methodOffset-1), ...
                        'RefAfter',refData(methodLast+1:methodLast+suffix(c)), ...
                        'SigBefore',sigData(methodOffset-prefix(c):methodOffset-1), ...
                        'SigAfter',sigData(methodLast+1:methodLast+suffix(c)));
                end
                run=measure(method,methodRef,methodSig,K,B,cfg,boundary);
                clear boundary
                run.Metric=section5_caf(methodRef,run.Residual,run.Info.OutputStart,interval,cfg);
                if methodFirst==cfg.startSample
                    run.RawMetric=results.Raw(intRow).Metric;
                else
                    run.RawMetric=section5_caf(methodRef,methodSig,1,interval,cfg);
                end
                if method=="ECA" || method=="ECA-FS", cache(key)=run; end
            end
            relativeError=NaN;
            if method=="ECA"
                referenceResidual=run.Residual;
            elseif method=="ECA-CA"
                denominator=norm(referenceResidual);
                if denominator>0
                    relativeError=100*norm(run.Residual-referenceResidual)/denominator;
                end
            end
            row=struct('CaseId',c,'Method',method,'N',N,'K',K,'B',B, ...
                'InBSweep',cases.InBSweep(c),'InKSweep',cases.InKSweep(c), ...
                'InNSweep',cases.InNSweep(c),'RelativeErrorPct',relativeError, ...
                'PTRdB',run.Metric.PTRdB,'MeanTimeSeconds',mean(run.Times), ...
                'TimingRepetitions',numel(run.Times), ...
                'PeakWorkingBytes',run.Info.PeakWorkingBytes, ...
                'PeakWorkingGB',run.Info.PeakWorkingBytes/1e9, ...
                'MinimumRcond',run.Info.MinimumRcond,'ResultReused',reused, ...
                'TargetMagnitude',run.Metric.TargetMagnitude, ...
                'TargetPower',run.Metric.TargetPower, ...
                'MeanSliceMagnitude',run.Metric.MeanSliceMagnitude, ...
                'ThresholdMagnitude',run.Metric.ThresholdMagnitude, ...
                'PTRDefinition',run.Metric.Definition, ...
                'PTRThresholdScale',run.Metric.ThresholdScale, ...
                'TargetIncludedInMean',run.Metric.TargetIncludedInMean, ...
                'RangeCellCount',run.Metric.RangeCellCount, ...
                'TotalRangeCellCount',run.Metric.TotalRangeCellCount, ...
                'MaskedRangeBinCount',numel(cfg.cafMaskedRangeBins), ...
                'DopplerMaskLowHz',cfg.cafDopplerMaskHz(1), ...
                'DopplerMaskHighHz',cfg.cafDopplerMaskHz(2), ...
                'CAFStart',interval(1),'CAFEnd',interval(2), ...
                'CAFDurationSeconds',run.Metric.DurationSeconds, ...
                'CAFReferenceBoundary',run.Metric.ReferenceBoundary, ...
                'CPIStartSample',methodFirst,'OutputLength',run.Info.OutputLength, ...
                'BaseCPIStartSample',cfg.startSample,'InputConvention',inputConvention, ...
                'SourceFrameIndex',sourceFrameIndex,'CAFReferenceStartSample',methodFirst, ...
                'DelayBoundary',run.Info.DelayBoundary,'CoefficientSolver',run.Info.CoefficientSolver, ...
                'InputLength',run.Info.InputLength, ...
                'InputStartSample',methodFirst+run.Info.InputStart-1, ...
                'PrefixSamples',run.Info.PrefixSamples,'SuffixSamples',run.Info.SuffixSamples, ...
                'BoundaryInputBytes',run.Info.BoundaryInputBytes, ...
                'DopplerBinSpacingHz',run.Metric.DopplerBinSpacingHz, ...
                'TargetDopplerHz',cfg.targetDopplerHz,'FFTLength',run.Info.FFTLength, ...
                'SegmentOrBlockCount',run.Info.Segments,'FSBatches',run.Info.FrequencyBatches);
            if isempty(rows), rows=row; else, rows(end+1)=row; end %#ok<AGROW>
            results.Slices{end+1}=run.Metric.Slice;
            detail=struct('CaseId',c,'Method',method,'CacheKey',key, ...
                'TimesSeconds',run.Times,'Memory',run.Info,'CAF',run.Metric, ...
                'RawCAF',run.RawMetric,'CPIStartSample',methodFirst, ...
                'InputStartSample',row.InputStartSample);
            if isempty(results.RunInfo)
                results.RunInfo=detail;
            else
                results.RunInfo(end+1)=detail;
            end
            fprintf('  %-7s PTR=%9.4f dB  mean=%9.4f s  memory=%.6f GB  reused=%d\n', ...
                method,row.PTRdB,row.MeanTimeSeconds,row.PeakWorkingGB,reused);
            clear run methodRef methodSig
        end
        results.Metrics=struct2table(rows);
        save_results(results);
    end
    results.Complete=true;
    save_results(results);
catch exception
    if ~isempty(rows), results.Metrics=struct2table(rows); end
    results.Failure=struct('Identifier',exception.identifier,'Message',exception.message);
    save_results(results);
    rethrow(exception);
end
if cfg.makePlots, plot_section5(fullfile(outDir,'results.mat')); end
fprintf('Saved %d configurations to %s\n',height(cases),outDir);
end

function vector = read_channel(file,name,first,count)
if ~isempty(whos('-file',file,'FirstSample'))
    metadata=load(file,'FirstSample');
    validateattributes(metadata.FirstSample,{'numeric'},{'scalar','integer','positive'});
    first=first-metadata.FirstSample+1;
end
description=whos('-file',file,name);
assert(isscalar(description) && any(description.size==1), ...
    'Expected a vector variable %s in %s.',name,file);
assert(first>=1 && first+count-1<=prod(description.size), ...
    'Requested CPI plus real ECA-FS boundary samples exceeds the recording.');
source=matfile(file);
if description.size(1)==1
    selection=substruct('.',name,'()',{1,first:first+count-1});
else
    selection=substruct('.',name,'()',{first:first+count-1,1});
end
vector=double(subsref(source,selection));
vector=complex(vector(:));
end

function run = measure(method,ref,sig,K,B,cfg,boundary)
for repetition=1:cfg.warmups
    warm=section5_cancel(method,ref,sig,K,B,cfg,false,boundary); %#ok<NASGU>
    clear warm
end
times=zeros(cfg.repetitions,1);
for repetition=1:cfg.repetitions
    timer=tic;
    residual=section5_cancel(method,ref,sig,K,B,cfg,false,boundary);
    times(repetition)=toc(timer);
    if repetition<cfg.repetitions, clear residual; end
end
% Instrumented call is outside all reported timings.
[checked,info]=section5_cancel(method,ref,sig,K,B,cfg,true,boundary);
assert(info.OutputStart==1 && numel(checked)==numel(ref), ...
    'Every method must return its complete retained N-sample interval.');
assert(norm(checked-residual)<=1e-12*max(norm(residual),1), ...
    'Memory instrumentation changed the cancellation result.');
run=struct('Residual',residual,'Times',times,'Info',info);
end

function save_results(results)
save(fullfile(results.OutputDirectory,'results.mat'),'results','-v7.3');
if ~isempty(results.Metrics)
    writetable(results.Metrics,fullfile(results.OutputDirectory,'metrics.csv'));
end
end
