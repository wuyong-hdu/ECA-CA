function metric = section5_apply_ptr(metric, cfg)
% Apply the common CAF mask and evaluate PTR over unmasked range cells.
% Shared by fresh CAF evaluation and re-evaluation of existing result files.
if ~isfield(metric,'UnmaskedSlice'), metric.UnmaskedSlice=metric.Slice; end
slice=metric.UnmaskedSlice;
assert(isrow(slice) && numel(slice)==cfg.rangeBins, ...
    'The saved slice must contain every configured range cell.');
assert(metric.DopplerHz==cfg.targetDopplerHz,'Saved target Doppler differs from the configuration.');
[slice,validRange]=section5_mask_caf(slice,1:cfg.rangeBins,metric.DopplerHz,cfg);
[ptr,meanMagnitude,~,threshold]=section5_ptr(slice,1,cfg.targetRangeBin, ...
    cfg.ptrThresholdScale,validRange);
metric.Slice=slice;
metric.ValidRangeMask=validRange;
metric.MeanRangeBins=find(validRange);
metric.CAFMask=struct('RangeBins',cfg.cafMaskedRangeBins, ...
    'DopplerHz',cfg.cafDopplerMaskHz,'FrequencyConvention',"fixed_hz_inclusive");
if isfield(metric,'LocalCAF')
    if ~isfield(metric,'UnmaskedLocalCAF'), metric.UnmaskedLocalCAF=metric.LocalCAF; end
    [metric.LocalCAF,metric.LocalValidMask]=section5_mask_caf( ...
        metric.UnmaskedLocalCAF,1:cfg.rangeBins,metric.LocalDopplerHz,cfg);
end
if isfield(metric,'RangeDopplerMap') && ~isempty(metric.RangeDopplerMap)
    map=metric.RangeDopplerMap;
    if ~isfield(map,'UnmaskedMagnitude'), map.UnmaskedMagnitude=map.Magnitude; end
    [map.Magnitude,map.ValidMask]=section5_mask_caf( ...
        map.UnmaskedMagnitude,map.RangeBins,map.DopplerHz,cfg);
    if isfield(map,'AllDopplerMagnitude')
        if ~isfield(map,'UnmaskedAllDopplerMagnitude')
            map.UnmaskedAllDopplerMagnitude=map.AllDopplerMagnitude;
        end
        [map.AllDopplerMagnitude,map.AllDopplerValidMask]=section5_mask_caf( ...
            map.UnmaskedAllDopplerMagnitude,map.RangeBins,map.AllDopplerHz,cfg);
    end
    map.CAFMask=metric.CAFMask;
    metric.RangeDopplerMap=map;
end
metric.PTRdB=ptr;
metric.TargetMagnitude=slice(cfg.targetRangeBin);
metric.TargetPower=metric.TargetMagnitude^2;
metric.MeanSliceMagnitude=meanMagnitude;
metric.ThresholdMagnitude=threshold;
metric.ThresholdScale=cfg.ptrThresholdScale;
metric.RangeCellCount=nnz(validRange);
metric.TotalRangeCellCount=numel(slice);
metric.TargetIncludedInMean=true;
metric.Definition="masked_slice_mean_magnitude_threshold_including_target";
end
