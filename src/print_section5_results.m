function print_section5_results(resultSource, caseIds)
% Print saved metrics in the original run_section5 console format.
% print_section5_results(resultFile)          % All saved configurations.
% print_section5_results(resultFile,[13 14])  % Original case numbers retained.
% Passing a file reloads retested timings; no algorithms or plots are run.
if isstruct(resultSource)
    r=resultSource;
else
    loaded=load(resultSource,'results'); r=loaded.results;
end
total=height(r.Cases);
if nargin<2,caseIds=1:total;end
if isempty(caseIds),return;end
validateattributes(caseIds,{'numeric'},{'vector','integer','positive','finite','<=',total});
caseIds=unique(caseIds,'stable');
methods=["ECA","ECA-B","ECA-FS","ECA-CA"];
for c=reshape(caseIds,1,[])
    fprintf('Case %d/%d: N=%d, K=%d, B=%d\n', ...
        c,total,r.Cases.N(c),r.Cases.K(c),r.Cases.B(c));
    for m=1:numel(methods)
        index=find(r.Metrics.CaseId==c & r.Metrics.Method==methods(m));
        if isempty(index),continue;end
        assert(isscalar(index),'Duplicate %s result in case %d.',methods(m),c);
        row=r.Metrics(index,:);
        fprintf('  %-7s PTR=%9.4f dB  mean=%9.4f s  memory=%.6f GB  reused=%d\n', ...
            methods(m),row.PTRdB,row.MeanTimeSeconds,row.PeakWorkingGB,row.ResultReused);
    end
    ca=r.Metrics(r.Metrics.CaseId==c & r.Metrics.Method=="ECA-CA",:);
    if ~isempty(ca)
        fprintf('  Relative residual error (ECA-CA/ECA): %.6e %%\n',ca.RelativeErrorPct);
    end
end
end
