function [ptrDB, meanMagnitude, cells, threshold] = section5_ptr(F, l, rd, scale, validRange)
% Peak-to-threshold ratio using unmasked cells of the target-Doppler slice.
% F contains magnitudes; validRange explicitly excludes artificial mask zeros.
% The target is included. Genuine zeros in valid cells still count in the mean.
if nargin<4, scale=4.5; end
if nargin<5, validRange=true(1,size(F,2)); end
validateattributes(F,{'double'},{'2d','real','nonnegative','finite','nonempty'});
validateattributes(scale,{'numeric'},{'scalar','real','finite','positive'});
validateattributes(l,{'numeric'},{'scalar','integer','positive'});
validateattributes(rd,{'numeric'},{'scalar','integer','positive'});
assert(l<=size(F,1) && rd<=size(F,2),'Target cell lies outside the CAF array.');
validateattributes(validRange,{'logical'},{'row','numel',size(F,2)});
assert(any(validRange) && validRange(rd),'section5:MaskedTarget', ...
    'The target and at least one averaging cell must lie outside the CAF mask.');
cells=F(l,validRange);
meanMagnitude=mean(cells);
threshold=scale*meanMagnitude;
ptrDB=20*log10(F(l,rd)/threshold);
end
