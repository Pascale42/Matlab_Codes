%
% function LoadCCG(FileBase, varargin)
%          [State Clust1 Clust2] = DefaultArgs(varargin,{[],[],[]});
%
% >>> Loads ccgs from mySpike2Spike.m and plot one, several or all the ccgs
% Clust 1&2 are cluster numbers, ex. [2 5] for shank 2 clu 5


function LoadCCG(FileBase, varargin);

[State Clust1 Clust2] = DefaultArgs(varargin,{[],[0 0],[0 0]});



%%% SPK PARAM %%%

if isempty(State)
    load([FileBase '.NeuronQuality.mat']);
else
    load([FileBase '.StateNeuronQuality' State '.mat']);
end
SpkPar=CatStruct(OutArgs);


%%% CCGs %%%

if isempty(State)
    load([FileBase '.s2s'],'-MAT');
else
    load([FileBase '.s2s.' State],'-MAT');
end

ElClu1=[];
ElClu2=[];

if Clust1 ~= 0
    for i=1:length(s2s.ElClu)
        if Clust1 == s2s.ElClu(i,:)
            ElClu1=i;
        end
    end
    if ElClu1 == 0
        error('Clust1 does not exist!')
        return
    end
    clear i;
end

if Clust2 ~= 0
    for i=1:length(s2s.ElClu)
        if Clust2 == s2s.ElClu(i,:)
            ElClu2=i;
        end
    end
    if isempty(ElClu2)
        error('Clust2 does not exist!')
        return
    end
    clear i;
end

if isempty(ElClu1) & isempty(ElClu2)

    figure(1)
    for i=1:length(s2s.ElClu)
        for j = 1:length(s2s.ElClu)
            subplot(131), bar(s2s.tbin,s2s.ccg(:,i,i)), colormap([0 0 0])
            axis tight; title(['Cluster ' num2str(s2s.ElClu(i,:)) '; ' num2str(SpkPar.FirRate(i)) ' Hz']);
            subplot(132), bar(s2s.tbin,s2s.ccg(:,i,j)), colormap([1 0 0])
            axis tight; title(['CCG ' num2str(s2s.ElClu(i,:)) ' vs ' num2str(s2s.ElClu(j,:))])
            subplot(133), bar(s2s.tbin,s2s.ccg(:,j,j)), colormap([0 0 0])
            axis tight; title(['Cluster ' num2str(s2s.ElClu(j,:)) '; ' num2str(SpkPar.FirRate(j)) ' Hz']);
            waitforbuttonpress
%             keydown = waitforbuttonpress
%             if (keydown == 1)
%                 disp('Key was pressed');
%             end
        end
    end
clear i;


elseif ~isempty(ElClu1) & isempty(ElClu2)
    figure(1)
    for i=1:length(s2s.ElClu)
        subplot(131), bar(s2s.tbin,s2s.ccg(:,ElClu1,ElClu1)), colormap([0 0 0])
        axis tight; title(['Cluster ' num2str(s2s.ElClu(ElClu1,:)) '; ' num2str(SpkPar.FirRate(ElClu1)) ' Hz']);
        subplot(132), bar(s2s.tbin,s2s.ccg(:,ElClu1,i)), colormap([1 0 0])
        axis tight; title(['CCG ' num2str(s2s.ElClu(ElClu1,:)) ' vs ' num2str(s2s.ElClu(i,:))])
        subplot(133), bar(s2s.tbin,s2s.ccg(:,i,i)), colormap([0 0 0])
        axis tight; title(['Cluster ' num2str(s2s.ElClu(i,:)) '; ' num2str(SpkPar.FirRate(i)) ' Hz']);
        waitforbuttonpress
    end

else
    figure(1)
    subplot(131), bar(s2s.tbin,s2s.ccg(:,ElClu1,ElClu1)), colormap([0 0 0])
    axis tight; title(['Cluster ' num2str(s2s.ElClu(ElClu1,:)) '; ' num2str(SpkPar.FirRate(ElClu1)) ' Hz']);
    subplot(132), bar(s2s.tbin,s2s.ccg(:,ElClu1,ElClu2)), colormap([1 0 0])
    axis tight; title(['CCG ' num2str(ElClu1) ' vs ' num2str(ElClu2)])
    subplot(133), bar(s2s.tbin,s2s.ccg(:,ElClu2,ElClu2)), colormap([0 0 0])
    axis tight; title(['Cluster ' num2str(s2s.ElClu(ElClu2,:)) '; ' num2str(SpkPar.FirRate(ElClu2)) ' Hz']);
end