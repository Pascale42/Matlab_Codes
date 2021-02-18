function PlotRaster(Res,Clu,varargin)
%function PlotRaster(Res,Clu,SampleRate)
[SampleRate] = DefaultArgs(varargin,{20000});

Colors = colormap;
Colors = repmat(Colors,10,1);
MyClu = unique(Clu);
%clf; hold on;
for c=1:length(MyClu)
    MyRes = Res(find(Clu==MyClu(c)));
    plot(repmat(MyRes(:)',2,1)/SampleRate, repmat(c+[0;.6],1,length(MyRes)), 'color', Colors(c,:));
    hold on
end    

set(gca, 'ytick', 1.45:length(MyClu)+.45);
set(gca, 'yticklabel', MyClu);
