function files = plot_section5_equivalence(resultSource)
% Numerical equivalence from the saved residual errors.
figureWidthCm = 17.5;
figureHeightCm = 5.8;
fontSizePt = 8;
fontName = 'Arial';
defaultYLimits = [1e-12 1e-10]; % Relative residual error in percent.

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
data=r.Metrics(r.Metrics.Method=="ECA-CA",:);
assert(~isempty(data),'No ECA-CA equivalence results were saved.');
errors=data.RelativeErrorPct;
assert(all(isnan(errors) | (isfinite(errors) & errors>=0)), ...
    'Relative errors must be nonnegative or undefined (NaN).');
positive=errors(isfinite(errors) & errors>0);
yLimits=defaultYLimits;
if ~isempty(positive)
    % Expand, never clip, when another saved run falls outside the defaults.
    yLimits(1)=min(yLimits(1),10^floor(log10(min(positive))));
    yLimits(2)=max(yLimits(2),10^ceil(log10(max(positive))));
end
valid=find(isfinite(errors));
maximumIndex=[];
if ~isempty(valid)
    [~,i]=max(errors(valid)); maximumIndex=valid(i);
end
groups={'InBSweep','InKSweep','InNSweep'};
parameters={'B','K','N'};
color=[0.20 0.20 0.20];
fig=figure('Visible','off','Units','centimeters', ...
    'Position',[2 2 figureWidthCm figureHeightCm], ...
    'PaperUnits','centimeters','PaperPositionMode','auto','Color','w');
cleanup=onCleanup(@() close_if_hidden(fig));
% Fixed equal-size axes keep horizontal labels and panel letters aligned.
leftMargin=0.09; rightMargin=0.02; panelGap=0.03;
panelWidth=(1-leftMargin-rightMargin-2*panelGap)/3;
panelBottom=0.29; panelHeight=0.66;
axesHandles=gobjects(1,3);
files=strings(0,1);
maxAnnotated=false;
for g=1:3
    values=sortrows(data(data.(groups{g}),:),parameters{g});
    x=values.(parameters{g});
    if g==3,x=x/1e6;end
    y=values.RelativeErrorPct;
    position=[leftMargin+(g-1)*(panelWidth+panelGap), ...
        panelBottom,panelWidth,panelHeight];
    ax=axes(fig,'Units','normalized','Position',position, ...
        'PositionConstraint','innerposition');
    axesHandles(g)=ax;
    hold(ax,'on');
    % A logarithmic axis with explicit doubling ticks has log2 spacing.
    set(ax,'XScale','log','YScale','log','YLim',yLimits, ...
        'FontName',fontName,'FontSize',fontSizePt, ...
        'Box','on','YMinorTick','off','YMinorGrid','off', ...
        'XMinorTick','off','XMinorGrid','off','XTickLabelRotation',0);
    yticks(ax,10.^(floor(log10(yLimits(1))):ceil(log10(yLimits(2)))));
    grid(ax,'on');
    visible=isfinite(y) & y>0;
    scatter(ax,x(visible),y(visible),24,color,'filled');
    if g==1
        ylabel(ax,'Relative residual error (%)', ...
            'FontName',fontName,'FontSize',fontSizePt);
    else
        ax.YTickLabel=[];
    end
    if g==3
        xTitle='{\it N} / 10^6';
    else
        xTitle=['{\it ' parameters{g} '}'];
    end
    xlabel(ax,xTitle,'Units','normalized','Position',[0.5 -0.22 0], ...
        'HorizontalAlignment','center','VerticalAlignment','baseline', ...
        'FontName',fontName,'FontSize',fontSizePt,'Interpreter','tex');
    if ~isempty(x)
        assert(all(isfinite(x) & x>0),'Logarithmic X values must be positive.');
        ticks=unique(x);
        xticks(ax,ticks);
        xticklabels(ax,compose('%g',ticks));
        xlim(ax,[min(x)/2^0.15 max(x)*2^0.15]);
    end
    annotation(fig,'textbox',[position(1),0.04,panelWidth,0.10], ...
        'String',sprintf('(%c)','a'+g-1),'LineStyle','none', ...
        'HorizontalAlignment','center','VerticalAlignment','middle', ...
        'FontName',fontName,'FontSize',fontSizePt,'Margin',0);
    if ~maxAnnotated && ~isempty(maximumIndex) && ...
            any(values.CaseId==data.CaseId(maximumIndex))
        pointIndex=find(values.CaseId==data.CaseId(maximumIndex),1);
        if y(pointIndex)>0
            % Convert the actual log-axis data point to figure coordinates.
            pointX=(log2(x(pointIndex))-log2(ax.XLim(1)))/diff(log2(ax.XLim));
            pointY=(log10(y(pointIndex))-log10(yLimits(1)))/diff(log10(yLimits));
            exponent=floor(log10(y(pointIndex)));
            mantissa=y(pointIndex)/10^exponent;
            annotation(fig,'textarrow', ...
                position(1)+position(3)*[0.95 pointX], ...
                position(2)+position(4)*[0.90 pointY], ...
                'String',sprintf('Max: %.2f\\times10^{%d}%%',mantissa,exponent), ...
                'Interpreter','tex','FontName',fontName,'FontSize',7.5, ...
                'HorizontalAlignment','right','VerticalAlignment','bottom', ...
                'HeadLength',4,'HeadWidth',4,'LineWidth',0.6);
        end
        maxAnnotated=true;
    end
    % Zero and undefined errors have no log coordinate; report explicitly.
    exceptional=strings(0,1);
    if any(y==0),exceptional(end+1)=sprintf('Zero error at %s: %s',parameters{g}, ...
            strtrim(sprintf('%g ',values.(parameters{g})(y==0))));end
    if any(isnan(y)),exceptional(end+1)=sprintf('Undefined at %s: %s',parameters{g}, ...
            strtrim(sprintf('%g ',values.(parameters{g})(isnan(y)))));end
    if ~isempty(exceptional)
        text(ax,0.02,0.04,strjoin(exceptional,newline),'Units','normalized', ...
            'FontName',fontName,'FontSize',7,'VerticalAlignment','bottom');
    end

end
linkaxes(axesHandles,'y');
base=fullfile(out,'fig3_equivalence');
section5_export_figure(fig,base,showFigures,true);
files(end+1,1)=string(base);
end

function close_if_hidden(fig)
if isgraphics(fig) && strcmp(fig.Visible,'off'),close(fig);end
end
