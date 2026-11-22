function layout = section5_fs_layout(N, B, cfg)
% Native FS segmentation needed to retain a complete N-sample output.
validateattributes(N,{'numeric'},{'scalar','integer','positive'});
validateattributes(B,{'numeric'},{'scalar','integer','positive'});
assert(mod(N,4*B)==0,'ECA-FS requires an integer step N/(4B).');
L=N/B; step=L/4;
validateattributes(cfg.fsOutputStart,{'numeric'},{'scalar','integer','positive'});
validateattributes(cfg.fsTailZeros,{'numeric'},{'scalar','integer','nonnegative','<',L});
assert(cfg.fsOutputStart+step-1<=L-cfg.fsTailZeros, ...
    'ECA-FS retained interval overlaps its zeroed tail.');
before=cfg.fsOutputStart-1;
after=L-step-before;
layout=struct('CPILength',N,'SegmentLength',L,'Step',step, ...
    'Segments',N/step,'PrefixSamples',before,'SuffixSamples',after, ...
    'InputLength',N+L-step,'InputStart',1-before,'OutputStart',1, ...
    'RetainedStart',cfg.fsOutputStart);
end
