function section5_export_figure(fig,base,showFigures,preservePageSize)
% Export a figure whose saved FIG opens visibly without a visibility override.
% showFigures controls whether the current session keeps its window open.
% preservePageSize exports the full physical page instead of tight cropping.
if nargin<3,showFigures=true;end
if nargin<4,preservePageSize=false;end
validateattributes(showFigures,{'logical'},{'scalar'});
validateattributes(preservePageSize,{'logical'},{'scalar'});
% MATLAB also consults the session default when opening a FIG without an
% explicit visibility option. Restore normal interactive opening when shown.
if showFigures,set(groot,'DefaultFigureVisible','on');end
if preservePageSize
    originalUnits=fig.Units;
    fig.Units='centimeters';
    pageSize=fig.Position(3:4);
    fig.Units=originalUnits;
    set(fig,'PaperUnits','centimeters','PaperSize',pageSize, ...
        'PaperPosition',[0 0 pageSize],'PaperPositionMode','manual');
    % Keep the whole page while embedding the selected TrueType fonts.
    % A white boundary prevents exportgraphics from tightly cropping it.
    pageFrame=annotation(fig,'rectangle',[0 0 1 1], ...
        'Color','w','LineWidth',0.01);
    exportgraphics(fig,[base '.pdf'],'ContentType','vector','BackgroundColor','white');
    delete(pageFrame);
    print(fig,[base '.png'],'-dpng','-r180');
else
    exportgraphics(fig,[base '.pdf'],'ContentType','vector');
    exportgraphics(fig,[base '.png'],'Resolution',180);
end
set(fig,'Visible','on');
drawnow;
savefig(fig,[base '.fig']);
if ~showFigures,close(fig);end
end
