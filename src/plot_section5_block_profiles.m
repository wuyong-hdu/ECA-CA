function files = plot_section5_block_profiles(resultSource, configurations)
% Fig. 6: four saved target-Doppler range profiles in a 2-by-2 layout.
% Default: N=1600000, K=400, B=[25 100 400 1600]. No algorithms are rerun.
% Each panel includes one gray before-suppression profile and four methods.
% Saved after-suppression PTR values are listed in two rows inside each panel.
% No threshold lines are plotted. Original FFT amplitudes are restored.
figureWidthCm=17.5;
figureHeightCm=11.0;
fontName='Arial';
fontSizePt=8;
if nargin<2
    configurations=[repmat([1600000 400],4,1),[25;100;400;1600]];
end
if isstruct(resultSource)
    r=resultSource;
else
    loaded=load(resultSource,'results'); r=loaded.results;
end
assert(r.Complete,'Only complete runs can generate comparison figures.');
validateattributes(configurations,{'numeric'}, ...
    {'size',[4 3],'integer','positive','finite'});
assert(size(unique(configurations,'rows'),1)==4, ...
    'Select four distinct [N K B] configurations.');
assert(all(configurations(:,1)==configurations(1,1)) && ...
    all(configurations(:,2)==configurations(1,2)), ...
    'The block-count comparison requires fixed N and K.');
[found,caseIds]=ismember(configurations,r.Cases{:,{'N','K','B'}},'rows');
assert(all(found),'Every selected configuration must exist in the saved results.');
showFigures=true;
if isfield(r.Config,'showFigures'),showFigures=r.Config.showFigures;end
out=fullfile(r.OutputDirectory,'figures');
if ~isfolder(out),mkdir(out);end
methods=["ECA","ECA-B","ECA-FS","ECA-CA"];
colors=[85 181 147;105 142 201;127 53 140;247 63 71]/255;
markers={'s','o','^','v'};
markerSizes=[4.5 4.0 4.2 3.8];
lineWidths=[1.5 0.9 0.9 0.8];
profiles=cell(4,5);
ptrs=zeros(4,4);
allFinite=[];
for p=1:4
    rawIndex=find([r.Raw.N]==configurations(p,1),1);
    assert(~isempty(rawIndex),'The before-suppression CAF was not saved.');
    profiles{p,1}=fft_magnitude_db(r.Raw(rawIndex).Metric,r.Config);
    for m=1:4
        runIndex=find([r.RunInfo.CaseId]==caseIds(p) & ...
            [r.RunInfo.Method]==methods(m),1);
        assert(~isempty(runIndex),'A requested method result was not saved.');
        profiles{p,m+1}=fft_magnitude_db(r.RunInfo(runIndex).CAF,r.Config);
        ptrs(p,m)=r.RunInfo(runIndex).CAF.PTRdB;
        assert(numel(profiles{p,m+1})==numel(profiles{p,1}), ...
            'All profiles must use the same range grid.');
    end
    for m=1:5
        values=profiles{p,m};
        allFinite=[allFinite;values(isfinite(values)).']; %#ok<AGROW>
    end
end
assert(~isempty(allFinite),'No finite CAF magnitudes were saved.');
yLimits=[20*floor(min(allFinite)/20),20*ceil(max(allFinite)/20)];
if yLimits(1)==yLimits(2),yLimits=yLimits+[-20 20];end
fig=figure('Visible','off','Units','centimeters', ...
    'Position',[2 2 figureWidthCm figureHeightCm], ...
    'PaperUnits','centimeters','PaperPositionMode','auto','Color','w');
cleanup=onCleanup(@() close_if_hidden(fig));
leftMargin=0.09; rightMargin=0.02; panelGap=0.06;
panelWidth=(1-leftMargin-rightMargin-panelGap)/2;
panelHeight=0.31;
rowBottoms=[0.60 0.13];
axesHandles=gobjects(1,4);
legendHandles=gobjects(1,5);
for p=1:4
    row=ceil(p/2); column=mod(p-1,2)+1;
    position=[leftMargin+(column-1)*(panelWidth+panelGap), ...
        rowBottoms(row),panelWidth,panelHeight];
    ax=axes(fig,'Units','normalized','Position',position, ...
        'PositionConstraint','innerposition');
    axesHandles(p)=ax;
    hold(ax,'on');
    set(ax,'FontName',fontName,'FontSize',fontSizePt, ...
        'Box','on','YLim',yLimits,'XTickLabelRotation',0);
    rangeBins=1:numel(profiles{p,1});
    beforeLine=plot(ax,rangeBins,profiles{p,1},'-', ...
        'Color',[0.55 0.55 0.55],'LineWidth',0.65, ...
        'DisplayName','Before clutter suppression');
    if p==1,legendHandles(1)=beforeLine;end
    % Sparse markers distinguish methods without obscuring the range profiles.
    markerIndices=unique([r.Config.targetRangeBin,20:20:(numel(rangeBins)-1)]);
    for m=1:4
        curve=plot(ax,rangeBins,profiles{p,m+1},'-','Color',colors(m,:), ...
            'Marker',markers{m},'MarkerIndices',markerIndices, ...
            'MarkerFaceColor',colors(m,:),'MarkerEdgeColor',colors(m,:), ...
            'MarkerSize',markerSizes(m),'LineWidth',lineWidths(m), ...
            'DisplayName',methods(m));
        if p==1,legendHandles(m+1)=curve;end
    end
    xline(ax,r.Config.targetRangeBin,':','Color',[0.25 0.25 0.25], ...
        'LineWidth',0.6,'HandleVisibility','off');
    xlim(ax,[0 numel(rangeBins)]);
    xticks(ax,0:20:numel(rangeBins));
    yticks(ax,(20*ceil(yLimits(1)/20)):20:yLimits(2));
    grid(ax,'on');
    % Use the clear lower band for a compact two-row PTR summary.
    text(ax,0.035,0.1225,'PTR (dB)','Units','normalized', ...
        'HorizontalAlignment','left','VerticalAlignment','middle', ...
        'FontName',fontName,'FontSize',fontSizePt,'Interpreter','none', ...
        'Color',[0.15 0.15 0.15],'BackgroundColor','w','Margin',0.5);
    ptrLefts=[0.35 0.69]; ptrRows=[0.18 0.065];
    for m=1:4
        ptrPosition=[ptrLefts(mod(m-1,2)+1),ptrRows(ceil(m/2))];
        text(ax,ptrPosition(1),ptrPosition(2), ...
            sprintf('%s: %.2f',methods(m),ptrs(p,m)), ...
            'Units','normalized','HorizontalAlignment','left', ...
            'VerticalAlignment','middle','FontName',fontName,'FontSize',fontSizePt, ...
            'Interpreter','none','Color',[0.15 0.15 0.15], ...
            'BackgroundColor','w','Margin',0.5, ...
            'Tag',sprintf('PTRValueCase%dMethod%d',p,m));
    end
    if column==1
        ylabel(ax,'CAF magnitude (dB)','FontName',fontName,'FontSize',fontSizePt);
    else
        ax.YTickLabel=[];
    end
    xlabel(ax,'Range bin','Units','normalized','Position',[0.5 -0.18 0], ...
        'HorizontalAlignment','center','VerticalAlignment','baseline', ...
        'FontName',fontName,'FontSize',fontSizePt);
    annotation(fig,'textbox',[position(1),position(2)-0.1075,panelWidth,0.05], ...
        'String',sprintf('(%c)','a'+p-1),'LineStyle','none', ...
        'HorizontalAlignment','center','VerticalAlignment','middle', ...
        'FontName',fontName,'FontSize',fontSizePt,'Margin',0);
end
linkaxes(axesHandles,'xy');
key=legend(axesHandles(1),legendHandles,["Before clutter suppression",methods], ...
    'Orientation','horizontal','NumColumns',5,'FontName',fontName, ...
    'FontSize',fontSizePt,'Box','on','AutoUpdate','off');
key.ItemTokenSize=[18 10];
key.Units='normalized';
drawnow;
keyPosition=key.Position;
keyPosition(1)=(leftMargin+1-rightMargin-keyPosition(3))/2;
keyPosition(2)=0.96-keyPosition(4)/2;
key.Position=keyPosition;
for p=1:4
    % Use the preserved ECA-CA peak to identify the common target range cell.
    add_target_label(fig,axesHandles(p), ...
        [r.Config.targetRangeBin,profiles{p,5}(r.Config.targetRangeBin)], ...
        [figureWidthCm figureHeightCm],fontName,fontSizePt);
end
base=fullfile(out,'fig6_block_range_profiles');
section5_export_figure(fig,base,showFigures,true);
files=string(base)+[".pdf";".png";".fig"];
fprintf('Saved block profiles at %.2f Hz; common magnitude limits [%g, %g] dB.\n', ...
    r.Config.targetDopplerHz,yLimits);
end

function values = fft_magnitude_db(caf,cfg)
assert(isfield(caf,'Slice') && isfield(caf,'DurationSeconds'), ...
    'The saved CAF slice and duration are required.');
assert(caf.DopplerHz==cfg.targetDopplerHz,'Target Doppler must match the saved slice.');
normalizer=round(caf.DurationSeconds*cfg.fs/cfg.cafAverageLength);
values=20*log10(abs(caf.Slice(:).')*normalizer);
values(~isfinite(values))=NaN; % Masked cells are not replaced by an artificial floor.
end

function close_if_hidden(fig)
if isgraphics(fig) && strcmp(fig.Visible,'off'),close(fig);end
end

function add_target_label(fig,ax,target,pageSize,fontName,fontSize)
position=ax.Position;
targetNormalized=(target-[ax.XLim(1),ax.YLim(1)])./[diff(ax.XLim),diff(ax.YLim)];
tail=position(1:2)+[0.40 0.89].*position(3:4);
tip=position(1:2)+targetNormalized.*position(3:4);
% Leave a small gap before the peak marker, matching Fig. 5 annotations.
direction=(tail-tip).*pageSize;
tip=tip+(0.08*direction/norm(direction))./pageSize;
annotation(fig,'textarrow',[tail(1) tip(1)],[tail(2) tip(2)], ...
    'String','Target','Color',[0.15 0.15 0.15],'LineWidth',0.6, ...
    'HeadLength',3,'HeadWidth',3,'FontName',fontName,'FontSize',fontSize, ...
    'Interpreter','none','HorizontalAlignment','left','VerticalAlignment','middle', ...
    'TextBackgroundColor','w','TextMargin',0.5);
end
