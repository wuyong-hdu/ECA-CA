function plot_section5(resultSource)
% Figures 3-8 of the paper, using saved metrics and CAF arrays only.
if isstruct(resultSource)
    r=resultSource;
else
    loaded=load(resultSource,'results'); r=loaded.results;
    r.OutputDirectory=fileparts(resultSource);
end
assert(r.Complete,'Only complete runs can generate comparison figures.');
plot_section5_equivalence(r);
plot_section5_ptr(r);
plot_section5_detectability(r,[1600000 400 25]);
plot_section5_block_profiles(r);
plot_section5_efficiency(r);
end
