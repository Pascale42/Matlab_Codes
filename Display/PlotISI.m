% 
% function PlotISI(FileBase)
% 

function PlotISI(FileBase)


[T,G,Map,Par]=LoadCluRes(FileBase);

[n, t] = ISIGrams(T, G);

for i=1:length(Map)
    bar(t,n(:,i))
    waitforbuttonpress
end