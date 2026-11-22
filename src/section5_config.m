function cfg = section5_config()
% Unified Section 5 settings. N is the complete output/evaluation CPI length.
% ECA-FS additionally reads real boundary samples from the same recording.
root = fileparts(mfilename('fullpath'));
cfg.dataFile = fullfile(root, '..', 'data', 'analysis_frame.mat');
cfg.referenceVariable = 'Ref';
cfg.surveillanceVariable = 'Sur';
cfg.startSample = 12 * 1600000 + 1; % Fixed recording segment for all N.
cfg.fs = 8e6;
cfg.targetRangeBin = 9;             % MATLAB bin 1 is zero delay.
cfg.targetDopplerHz = -55;
cfg.rangeBins = 100;
cfg.ptrThresholdScale = 4.5;       % Threshold = scale * mean of unmasked range magnitudes.
cfg.cafMaskedRangeBins = 1;        % TR=1: first range cell is excluded, including from the mean.
cfg.cafDopplerMaskHz = [-5 5];     % Fixed physical-frequency interval for every CPI.
cfg.cafDopplerHalfWidth = 4;       % Auxiliary saved neighborhood; PTR uses its center only.
cfg.cafAverageLength = 1000;
cfg.collectCAFMap = true;         % Save a 2D CAF for every evaluated output.
cfg.cafMapRangeLimits = [0 100];   % Display limits; data retain bins 1:100.
cfg.cafMapDopplerLimitsHz = [-150 150];
cfg.cafMapFFTColorLimitsDB = [];  % Raw FFT dB; empty = independent automatic scale per panel.
cfg.warmups = 1;
cfg.repetitions = 50;
cfg.fsTailZeros = 90;
cfg.fsOutputStart = 161;           % One-based index within every FS segment.
cfg.fsBatchTarget = 1000;
cfg.fsExtendedBatch = 1050;
cfg.fsLambda = 1e-10;
cfg.fsInputConvention = 'legacy_sliding'; % Original frame*extended-length acquisition.
cfg.ecabImplementation = 'legacy_batch';  % Batch-circular delays, backslash, serial workspace.
cfg.makePlots = true;
cfg.showFigures = true;          % Display generated plots and keep their windows open.
cfg.outputRoot = fullfile(root, '..', 'runs');
cfg.label = 'section5';
% Each row is [N K B]. The common case is evaluated once.
bCases = [repmat([1600000 400], 8, 1), [25;50;100;200;400;800;1600;3200]];
kCases = [repmat(1600000,4,1), [100;200;400;800], repmat(200,4,1)];
nCases = [[200000;400000;800000;1600000], repmat([400 200],4,1)];
cases = unique([bCases;kCases;nCases], 'rows', 'stable');
cfg.cases = array2table(cases, 'VariableNames', {'N','K','B'});
cfg.cases.InBSweep = ismember(cases,bCases,'rows');
cfg.cases.InKSweep = ismember(cases,kCases,'rows');
cfg.cases.InNSweep = ismember(cases,nCases,'rows');
end
