function raster(in)


if size(in,1) > size(in,2)
    in=in';
end
axis([0 max(in)+1 -1 2])
plot([in;in],[ones(size(in));zeros(size(in))],'k-')
set(gca,'TickDir','out') % draw the tick marks on the outside
set(gca,'YTick', []) % don't draw y-axis ticks
set(gca,'PlotBoxAspectRatio',[1 0.05 1]) % short and wide
set(gca,'Color',get(gcf,'Color')) % match figure background
set(gca,'YColor',get(gcf,'Color')) % hide the y axis
box off
