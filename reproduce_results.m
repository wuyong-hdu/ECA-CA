function results = reproduce_results(showFigures)
% Display the published results and recreate Figures 3-8 without raw data.
if nargin<1,showFigures=true;end
validateattributes(showFigures,{'logical'},{'scalar'});
root=fileparts(mfilename('fullpath'));
previousPath=path;
cleanup=onCleanup(@() path(previousPath));
addpath(fullfile(root,'src'));
loaded=load(fullfile(root,'results','paper_results.mat'),'results');
results=loaded.results;
results.OutputDirectory=fullfile(root,'output');
results.Config.showFigures=showFigures;
if ~isfolder(results.OutputDirectory),mkdir(results.OutputDirectory);end
report=evalc('print_section5_results(results);');
fprintf('%s',report);
fid=fopen(fullfile(results.OutputDirectory,'summary.txt'),'w');
assert(fid~=-1,'Cannot write summary.txt.');
fileCleanup=onCleanup(@() fclose(fid));
fprintf(fid,'%s',report);
clear fileCleanup
writetable(results.Metrics,fullfile(results.OutputDirectory,'metrics.csv'));
plot_section5(results);
fprintf('Figures and tables: %s\n',results.OutputDirectory);
end
