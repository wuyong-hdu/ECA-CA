function results = run_experiments(dataFile, cfg)
% Recompute the 14 configurations using a MAT file containing Ref and Sur.
root=fileparts(mfilename('fullpath'));
previousPath=path;
pathCleanup=onCleanup(@() path(previousPath));
addpath(fullfile(root,'src'));
if nargin<2,cfg=section5_config();end
if nargin<1,dataFile=cfg.dataFile;end
assert(isfile(dataFile),'Provide the acquisition MAT file containing Ref and Sur.');
cfg.dataFile=dataFile;
previousThreads=maxNumCompThreads;
threadCleanup=onCleanup(@() maxNumCompThreads(previousThreads));
maxNumCompThreads(8);
results=run_section5(cfg);
end
