function files = plot_section5_ptr(resultSource)
% Replot saved detectability metrics without running algorithms or timing.
figureWidthCm=17.5;
figureHeightCm=5.8;
fontName='Arial';
fontSizePt=8;
if isstruct(resultSource)
    r=resultSource;
else
    loaded=load(resultSource,'results'); r=loaded.results;
end
assert(r.Complete,'Only complete runs can generate comparison figures.');
showFigures=true;
if isfield(r.Config,'showFigures'),showFigures=r.Config.showFigures;end
out=fullfile(r.OutputDirectory,'figures');
if ~isfolder(out),mkdir(out);end

metricField='PTRdB'; metricLabel='PTR (dB)'; metricFigure='fig4_ptr';
methods=["ECA","ECA-B","ECA-FS","ECA-CA"];
colors=[85 181 147;105 142 201;127 53 140;247 63 71]/255;
markers={'s','o','^','v'};
% Distinct marker sizes and line widths expose the coincident ECA/CA curves.
markerSizes=[5.5 4.8 5.0 4.4];
lineWidths=[2.0 1.2 1.2 1.0];
groups={'InBSweep','InKSweep','InNSweep'};
parameters={'B','K','N'};
metricValues=r.Metrics.(metricField);
finiteValues=metricValues(isfinite(metricValues));
assert(~isempty(finiteValues),'No finite detectability values were saved.');
yLimits=[10*floor(min(finiteValues)/10),10*ceil(max(finiteValues)/10)];
if yLimits(1)==yLimits(2),yLimits=yLimits+[-10 10];end
fig=figure('Visible','off','Units','centimeters', ...
    'Position',[2 2 figureWidthCm figureHeightCm], ...
    'PaperUnits','centimeters','PaperPositionMode','auto','Color','w');
cleanup=onCleanup(@() close_if_hidden(fig));
leftMargin=0.09; rightMargin=0.02; panelGap=0.03;
panelWidth=(1-leftMargin-rightMargin-2*panelGap)/3;
panelBottom=0.29; panelHeight=0.59;
axesHandles=gobjects(1,3);
legendHandles=gobjects(1,4);
for g=1:3
    position=[leftMargin+(g-1)*(panelWidth+panelGap), ...
        panelBottom,panelWidth,panelHeight];
    ax=axes(fig,'Units','normalized','Position',position, ...
        'PositionConstraint','innerposition');
    axesHandles(g)=ax;
    hold(ax,'on');
    set(ax,'XScale','log','YLim',yLimits,'FontName',fontName, ...
        'FontSize',fontSizePt,'Box','on','XTickLabelRotation',0, ...
        'XMinorTick','off','XMinorGrid','off','YMinorTick','off');
    yticks(ax,yLimits(1):10:yLimits(2));
    grid(ax,'on');
    panelRows=r.Metrics.(groups{g});
    panelValues=r.Metrics(panelRows,:);
    allX=unique(panelValues.(parameters{g}));
    if g==3,allX=allX/1e6;end
    assert(~isempty(allX) && all(isfinite(allX) & allX>0), ...
        'Logarithmic X values must be positive.');
    for m=1:4
        values=sortrows(panelValues(panelValues.Method==methods(m),:),parameters{g});
        x=values.(parameters{g});
        if g==3,x=x/1e6;end
        curve=plot(ax,x,values.(metricField),'-','Color',colors(m,:), ...
            'Marker',markers{m},'MarkerFaceColor',colors(m,:), ...
            'MarkerEdgeColor',colors(m,:),'MarkerSize',markerSizes(m), ...
            'LineWidth',lineWidths(m),'DisplayName',methods(m));
        if g==1,legendHandles(m)=curve;end
    end
    xticks(ax,allX);
    xticklabels(ax,compose('%g',allX));
    xlim(ax,[min(allX)/2^0.15 max(allX)*2^0.15]);
    if g==1
        ylabel(ax,metricLabel,'FontName',fontName,'FontSize',fontSizePt);
    else
        ax.YTickLabel=[];
    end
    xTitle=['{\it ' parameters{g} '}'];
    if g==3,xTitle='{\it N} / 10^6';end
    xlabel(ax,xTitle,'Units','normalized','Position',[0.5 -0.246 0], ...
        'HorizontalAlignment','center','VerticalAlignment','baseline', ...
        'FontName',fontName,'FontSize',fontSizePt,'Interpreter','tex');
    annotation(fig,'textbox',[position(1),0.04,panelWidth,0.10], ...
        'String',sprintf('(%c)','a'+g-1),'LineStyle','none', ...
        'HorizontalAlignment','center','VerticalAlignment','middle', ...
        'FontName',fontName,'FontSize',fontSizePt,'Margin',0);
end
linkaxes(axesHandles,'y');
key=legend(axesHandles(1),legendHandles,methods,'Orientation','horizontal', ...
    'NumColumns',4,'FontName',fontName,'FontSize',fontSizePt, ...
    'Box','on','AutoUpdate','off');
key.ItemTokenSize=[18 10];
key.Units='normalized';
drawnow;
keyPosition=key.Position;
keyPosition(1)=(leftMargin+1-rightMargin-keyPosition(3))/2;
keyPosition(2)=0.95-keyPosition(4)/2;
key.Position=keyPosition;
base=fullfile(out,metricFigure);
section5_export_figure(fig,base,showFigures,true);
files=string(base)+[".pdf";".png";".fig"];
end

function close_if_hidden(fig)
if isgraphics(fig) && strcmp(fig.Visible,'off'),close(fig);end
end
