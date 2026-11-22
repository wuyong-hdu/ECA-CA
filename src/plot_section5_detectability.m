function files = plot_section5_detectability(resultSource, configuration)
% Fig. 5: before/after range-Doppler maps and target slices, four method rows.
% Uses one shared before-suppression CAF and the saved after-suppression CAFs.
% Default [N K B]=[1600000 400 25]. No signal processing or timing is rerun.
figureWidthCm=17.5;
figureHeightCm=20.5;
fontName='Arial';
fontSizePt=8;
rangeTicks=0:20:100;
if nargin<2,configuration=[1600000 400 25];end
validateattributes(configuration,{'numeric'}, ...
    {'size',[1 3],'integer','positive','finite'});
if isstruct(resultSource)
    r=resultSource;
else
    loaded=load(resultSource,'results'); r=loaded.results;
end
assert(r.Complete,'Only complete runs can generate comparison figures.');
[found,caseId]=ismember(configuration,r.Cases{:,{'N','K','B'}},'rows');
assert(found,'The selected [N K B] must exist in the saved results.');
showFigures=true;
if isfield(r.Config,'showFigures'),showFigures=r.Config.showFigures;end
out=fullfile(r.OutputDirectory,'figures');
if ~isfolder(out),mkdir(out);end
methods=["ECA","ECA-B","ECA-FS","ECA-CA"];
colors=[85 181 147;105 142 201;127 53 140;247 63 71]/255;
markers={'s','o','^','v'};
rawIndex=find([r.Raw.N]==configuration(1),1);
assert(~isempty(rawIndex),'The before-suppression CAF was not saved.');
raw=r.Raw(rawIndex).Metric;
rawMap=map_display(raw.RangeDopplerMap,r.Config);
rawSlice=slice_display(raw,r.Config);
afterMaps=cell(1,4);
afterSlices=cell(1,4);
thresholds=zeros(1,4);
ptrs=zeros(1,4);
allProfileValues=rawSlice(isfinite(rawSlice));
for m=1:4
    index=find([r.RunInfo.CaseId]==caseId & [r.RunInfo.Method]==methods(m),1);
    assert(~isempty(index),'A requested method result was not saved.');
    caf=r.RunInfo(index).CAF;
    afterMaps{m}=map_display(caf.RangeDopplerMap,r.Config);
    afterSlices{m}=slice_display(caf,r.Config);
    assert(numel(afterSlices{m})==numel(rawSlice),'CAF range grids must match.');
    normalizer=round(caf.DurationSeconds*r.Config.fs/r.Config.cafAverageLength);
    thresholds(m)=20*log10(caf.ThresholdMagnitude*normalizer);
    assert(isfinite(thresholds(m)),'The saved threshold must be finite and positive.');
    ptrs(m)=caf.PTRdB;
    values=afterSlices{m};
    allProfileValues=[allProfileValues,values(isfinite(values)),thresholds(m)]; %#ok<AGROW>
end
profileLimits=[20*floor(min(allProfileValues)/20),20*ceil(max(allProfileValues)/20)];
if profileLimits(1)==profileLimits(2),profileLimits=profileLimits+[-20 20];end
fig=figure('Visible','off','Units','centimeters', ...
    'Position',[2 2 figureWidthCm figureHeightCm], ...
    'PaperUnits','centimeters','PaperPositionMode','auto','Color','w');
cleanup=onCleanup(@() close_if_hidden(fig));
% Physical positions keep the twelve panels legible at full paper width.
mapLefts=[1.5 6.6]; mapWidth=3.6;
profileLeft=12.0; profileWidth=5.15;
panelHeight=3.3; rowBottoms=[15.8 11.05 6.3 1.55];
pageSize=[figureWidthCm figureHeightCm];
header_centers=[mapLefts+mapWidth/2,profileLeft+profileWidth/2];
headerGapCm=0.2;
headerBottomCm=rowBottoms(1)+panelHeight+headerGapCm;
headers={'Before clutter suppression','After clutter suppression', ...
    sprintf('Range profile at %.0f Hz',r.Config.targetDopplerHz)};
for col=1:3
    add_text(fig,[header_centers(col)-2.65,headerBottomCm,5.3,0.4],headers{col}, ...
        pageSize,fontName,fontSizePt,'normal');
end
mapAxes=gobjects(4,2);
profileAxes=gobjects(1,4);
for m=1:4
    bottom=rowBottoms(m);
    for col=1:2
        position=[mapLefts(col),bottom,mapWidth,panelHeight];
        ax=axes(fig,'Units','centimeters','Position',position, ...
            'PositionConstraint','innerposition');
        mapAxes(m,col)=ax;
        data=rawMap;
        if col==2,data=afterMaps{m};end
        imagesc(ax,data.RangeBins,data.Frequencies,data.DB,data.ColorLimits);
        axis(ax,'xy');
        xlim(ax,[0 100]); ylim(ax,[-150 150]);
        colormap(ax,parula(256)); hold(ax,'on');
        if col==2
            plot(ax,r.Config.targetRangeBin,r.Config.targetDopplerHz,'wo', ...
                'MarkerSize',4.5,'LineWidth',0.9,'HandleVisibility','off');
        end
        set(ax,'FontName',fontName,'FontSize',fontSizePt,'Box','on', ...
            'TickDir','out','Layer','top','XTickLabelRotation',0);
        xticks(ax,rangeTicks); yticks(ax,-150:50:150);
        if col==1
            ylabel(ax,'Doppler (Hz)','FontName',fontName,'FontSize',fontSizePt);
            text(ax,-0.375,0.5,methods(m),'Units','normalized', ...
                'Rotation',90,'HorizontalAlignment','center', ...
                'VerticalAlignment','middle','FontName',fontName, ...
                'FontSize',fontSizePt,'FontWeight','bold','Clipping','off');
        else
            ax.YTickLabel=[];
        end
        xlabel(ax,'Range bin','Units','normalized','Position',[0.5 -0.197 0], ...
            'HorizontalAlignment','center','VerticalAlignment','baseline', ...
            'FontName',fontName,'FontSize',fontSizePt);
        cb=colorbar(ax);
        set(cb,'Units','centimeters','Position', ...
            [position(1)+mapWidth+0.10,bottom,0.16,panelHeight], ...
            'FontName',fontName,'FontSize',fontSizePt);
        ax.Position=position;
        if col==2
            add_target_label(fig,ax,[0.38 0.15], ...
                [r.Config.targetRangeBin,r.Config.targetDopplerHz], ...
                pageSize,fontName,fontSizePt);
        end
        add_panel_letter(fig,position,(m-1)*3+col,pageSize,fontName,fontSizePt);
    end
    position=[profileLeft,bottom,profileWidth,panelHeight];
    ax=axes(fig,'Units','centimeters','Position',position, ...
        'PositionConstraint','innerposition');
    profileAxes(m)=ax;
    hold(ax,'on');
    rangeBins=1:numel(rawSlice);
    plot(ax,rangeBins,rawSlice,'-','Color',[0.55 0.55 0.55], ...
        'LineWidth',0.65,'Tag','SharedBeforeProfile');
    markerIndices=unique([r.Config.targetRangeBin,20:20:(numel(rangeBins)-1)]);
    plot(ax,rangeBins,afterSlices{m},'-','Color',colors(m,:), ...
        'LineWidth',1.0,'Marker',markers{m},'MarkerIndices',markerIndices, ...
        'MarkerSize',3.5,'MarkerFaceColor',colors(m,:), ...
        'MarkerEdgeColor',colors(m,:),'Tag','AfterProfile');
    yline(ax,thresholds(m),'--','Color',colors(m,:), ...
        'LineWidth',0.9,'Tag','PTRThresholdAfter');
    xline(ax,r.Config.targetRangeBin,':','Color',[0.25 0.25 0.25], ...
        'LineWidth',0.6,'HandleVisibility','off');
    set(ax,'FontName',fontName,'FontSize',fontSizePt,'Box','on', ...
        'XLim',[0 100],'YLim',profileLimits,'XTickLabelRotation',0);
    xticks(ax,rangeTicks); yticks(ax,profileLimits(1):20:profileLimits(2));
    grid(ax,'on');
    xlabel(ax,'Range bin','Units','normalized','Position',[0.5 -0.197 0], ...
        'HorizontalAlignment','center','VerticalAlignment','baseline', ...
        'FontName',fontName,'FontSize',fontSizePt);
    ylabel(ax,'CAF magnitude (dB)','FontName',fontName,'FontSize',fontSizePt);
    add_target_label(fig,ax,[0.42 0.89], ...
        [r.Config.targetRangeBin,afterSlices{m}(r.Config.targetRangeBin)], ...
        pageSize,fontName,fontSizePt);
    text(ax,0.96,0.10,sprintf('PTR = %.2f dB',ptrs(m)), ...
        'Units','normalized','HorizontalAlignment','right', ...
        'VerticalAlignment','middle','FontName',fontName,'FontSize',fontSizePt, ...
        'Color',[0.15 0.15 0.15],'BackgroundColor','w','Margin',1, ...
        'Interpreter','none','Tag','PTRValueAfter');
    add_panel_letter(fig,position,m*3,pageSize,fontName,fontSizePt);
end
linkaxes(profileAxes,'xy');
% The paper caption explains before/after curves and the threshold line.
drawnow;
% Confirm the displayed before maps and profiles are identical across rows.
for m=2:4
    firstImage=findobj(mapAxes(1,1),'Type','image');
    thisImage=findobj(mapAxes(m,1),'Type','image');
    assert(isequaln(firstImage.CData,thisImage.CData) && ...
        isequal(mapAxes(1,1).CLim,mapAxes(m,1).CLim),'Before maps must match.');
end
for m=1:4
    margin=afterSlices{m}(r.Config.targetRangeBin)-thresholds(m);
    assert(abs(margin-ptrs(m))<1e-9,'The displayed profile must preserve the saved PTR.');
    fprintf('%s: displayed target-to-threshold margin %.6f dB.\n',methods(m),margin);
end
base=fullfile(out,'fig5_detectability');
section5_export_figure(fig,base,showFigures,true);
files=string(base)+[".pdf";".png";".fig"];
end

function data = map_display(map,cfg)
assert(~isempty(map),'Saved range-Doppler maps are required.');
normalizer=round(cfg.fs/(cfg.cafAverageLength*map.DopplerBinSpacingHz));
if isfield(map,'AllDopplerMagnitude')
    magnitude=map.AllDopplerMagnitude; frequencies=map.AllDopplerHz;
else
    magnitude=map.Magnitude; frequencies=map.DopplerHz;
end
db=20*log10(double(magnitude)*normalizer);
finite=db(isfinite(db));
assert(~isempty(finite),'The saved map has no finite magnitudes.');
limits=[min(finite),max(finite)];
if limits(1)==limits(2),limits=limits+[-1 1];end
if isfield(cfg,'cafMapFFTColorLimitsDB') && ~isempty(cfg.cafMapFFTColorLimitsDB)
    limits=cfg.cafMapFFTColorLimitsDB;
end
data=struct('DB',db,'Frequencies',frequencies, ...
    'RangeBins',map.RangeBins,'ColorLimits',limits);
end

function values = slice_display(caf,cfg)
assert(caf.DopplerHz==cfg.targetDopplerHz,'The saved target Doppler must match.');
normalizer=round(caf.DurationSeconds*cfg.fs/cfg.cafAverageLength);
values=20*log10(abs(caf.Slice(:).')*normalizer);
values(~isfinite(values))=NaN;
end

function add_panel_letter(fig,position,index,pageSize,fontName,fontSize)
add_text(fig,[position(1),position(2)-1.23,position(3),0.4], ...
    sprintf('(%c)','a'+index-1),pageSize,fontName,fontSize,'normal');
end

function add_target_label(fig,ax,labelPosition,target,pageSize,fontName,fontSize)
position=ax.Position;
targetNormalized=(target-[ax.XLim(1),ax.YLim(1)])./[diff(ax.XLim),diff(ax.YLim)];
tail=position(1:2)+labelPosition.*position(3:4);
tip=position(1:2)+targetNormalized.*position(3:4);
% End just short of the marker so the target response stays visible.
direction=tail-tip;
tip=tip+0.08*direction/norm(direction);
annotation(fig,'textarrow',[tail(1) tip(1)]/pageSize(1), ...
    [tail(2) tip(2)]/pageSize(2),'String','Target', ...
    'Color',[0.15 0.15 0.15],'LineWidth',0.6,'HeadLength',3,'HeadWidth',3, ...
    'FontName',fontName,'FontSize',fontSize,'Interpreter','none', ...
    'HorizontalAlignment','left','VerticalAlignment','middle', ...
    'TextBackgroundColor','w','TextMargin',0.5);
end

function add_text(fig,position,label,pageSize,fontName,fontSize,fontWeight)
normalized=position./[pageSize pageSize];
annotation(fig,'textbox',normalized,'String',label,'LineStyle','none', ...
    'HorizontalAlignment','center','VerticalAlignment','middle', ...
    'FontName',fontName,'FontSize',fontSize,'FontWeight',fontWeight, ...
    'Margin',0,'Interpreter','none');
end

function close_if_hidden(fig)
if isgraphics(fig) && strcmp(fig.Visible,'off'),close(fig);end
end
