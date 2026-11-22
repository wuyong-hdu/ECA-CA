function intervals = section5_input_intervals(N,B,cfg)
% Actual sample intervals for the native sliding FS implementation.
layout=section5_fs_layout(N,B,cfg);
convention="cpi_boundary";
if isfield(cfg,'fsInputConvention'),convention=string(cfg.fsInputConvention);end
switch convention
    case "legacy_sliding"
        % Express the selected base start as the corresponding N-sample frame.
        % This keeps startSample fixed across the N sweep and matches frame=12
        % in the original N=1600000 example.
        assert(mod(cfg.startSample-1,N)==0, ...
            'Legacy sliding acquisition requires startSample-1 to be divisible by N.');
        frameIndex=(cfg.startSample-1)/N;
        inputFirst=frameIndex*layout.InputLength+1;
        outputFirst=inputFirst+layout.PrefixSamples;
    case "cpi_boundary"
        frameIndex=NaN;
        outputFirst=cfg.startSample;
        inputFirst=outputFirst-layout.PrefixSamples;
    otherwise
        error('Unknown FS input convention: %s',convention);
end
intervals=struct('Convention',convention,'FrameIndex',frameIndex, ...
    'BaseStart',cfg.startSample,'BaseEnd',cfg.startSample+N-1, ...
    'FSInputStart',inputFirst,'FSInputEnd',inputFirst+layout.InputLength-1, ...
    'FSOutputStart',outputFirst,'FSOutputEnd',outputFirst+N-1);
end
