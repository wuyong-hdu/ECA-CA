function [masked, valid] = section5_mask_caf(magnitude, rangeBins, dopplerHz, cfg)
% Common CAF postprocessing mask; zeros inside the mask are not observations.
% Keep genuine zero-valued samples outside the mask in all averaging counts.
validateattributes(magnitude,{'double'},{'2d','real','nonnegative','finite'});
assert(size(magnitude,1)==numel(dopplerHz) && size(magnitude,2)==numel(rangeBins));
validateattributes(cfg.cafMaskedRangeBins,{'numeric'},{'vector','integer','positive','<=',cfg.rangeBins});
validateattributes(cfg.cafDopplerMaskHz,{'numeric'},{'vector','numel',2,'real','finite'});
limits=cfg.cafDopplerMaskHz;
assert(limits(1)<=limits(2));
rangeKeep=~ismember(reshape(rangeBins,1,[]),cfg.cafMaskedRangeBins);
frequencyKeep=~(dopplerHz(:)>=limits(1) & dopplerHz(:)<=limits(2));
valid=frequencyKeep & rangeKeep;
masked=magnitude;
masked(~valid)=0;
end
