%
% A simple matlab function to load from the .clu and .res files the vectors G and T
% from many electrodes
%
% [T,G,Map,Par]=LoadCluRes(FileBase, | ElGpsToLoad,SaveFile, CluLocOnly, ClusToLoad)
% DefaultArgs(varargin,{[1:Par.nElecGps],0,[],[]});
%
% ElGpsToLoad; loaded from .xml by default, otherwise precise [1 2 ..]
% T is in samples
% CluLocOnly ; allows to specify the channels from which the clusters are
% loaded; you must specify ElGpsToLoad; it will also run ClusterLocation accordingly
% ClusToLoad allows to specifiy the vector of El numbers in ElGpsToLoad and
% matching length vector of Clu numbers (original indexing) such that only
% those are loaded
% alternatively ClusToLoad can specify ElClu pairs , and then ElGspToLoad
% can be just unique list of electrodes to load where these ElClu pairs are
% T is in samples, G goes from 1 to total number of clusters (excludes
% noise and artifacts )
% Map is a matrix displaying the correspondence between new cluster numbers (first column) and inital
% shank # (second column) and cluster # (third column)
%
% Pascale production, Anton just helped :))

function [T,G,Map,Par]=LoadCluRes(FileBase,varargin)

Par = LoadPar([FileBase '.xml']);

[ElGpsToLoad,SaveFile, CluLocOnly, ClusToLoad] = DefaultArgs(varargin,{1:Par.nElecGps,1,[],[]});

G=[];
T=[];
Map=[];
maxG=0;


%%% Loop over x=ElGpsToLoad from LoadPar

for x=ElGpsToLoad(:)'            % for each ElGp, load clu and res
    
    if ~FileExists([FileBase '.clu.' num2str(x)])
        continue
    end
    g=LoadClu([FileBase '.clu.' num2str(x)]);
%     fid = fopen([FileBase '.res.' num2str(x)],'r');
%     t = fscanf(fid,'%d');        %  or 
    t=load([FileBase '.res.' num2str(x)]);
    
    %%% Removes clusters artifact and noise clusters (0 & 1)
    indx=find(g>1);
    if isempty(indx);
        continue
    end
    g=g(indx);
    t=t(indx);
    
    %%% Creates vector of initial g and renames cluster # since 1 to n
    [ugini,~,g]=unique(g); % ugini: vector of initial unique values of g
    g=maxG+g;
    ug=unique(g);
    
    %%% Concatenates all the g and t
    G=[G;g];
    maxG=max(G);
    T=[T;t];
    
    %%% Creates a "map" matrix
    shk=zeros(length(ug),1)+x; % shank #
    map=[ug,shk,ugini];
    Map=[Map;map];
end

%%% Sort the spikes not to have surprizes later -A
[T, si] = sort(T);
G=G(si);

%%% Now more fancy : if one specifies [El Clu] pairs to load
if ~isempty(ClusToLoad)
    if length(ClusToLoad)~=length(ElGpsToLoad) && size(ClusToLoad,2)~=2
        warning('length(ClusToLoad)~=length(ElGpsToLoad)');
        return;
    end
    if size(ClusToLoad,2)~=2
        myi = find(ismember(Map(:,2:3),[ElGpsToLoad(:), ClusToLoad(:)],'rows'));
    else
        myi = find(ismember(Map(:,2:3),ClusToLoad,'rows'));
    end
    gi = ismember(G, myi);
    G = G(gi);
    T = T(gi);
    [~ , ~, G] = unique(G);
    Map = [1:length(myi)' Map(myi,2:3)]; %#ok<BDSCA>
    
end

%%% Even more fancy, load from specific channels according to ClusterLocation  # 30 avril 2015 #
if ~isempty(CluLocOnly)
    cluloc=ClusterLocation(FileBase, ElGpsToLoad,0);
    cluloc=[Map(:,1) cluloc];
    
    GoodClu=[];
    for n=1:size(cluloc,1)
        if ismember(cluloc(n,4),CluLocOnly)
            GoodClu=[GoodClu;cluloc(n,:)];
        end
    end
    
    Map=GoodClu(:,1:4);
    
    Indx=[];
    for n=1:size(GoodClu,1)
       indx=find(G == GoodClu(n,1));
        Indx=[Indx; indx];
    end
    Indx=sort(Indx);
    G=G(Indx);
    T=T(Indx);
    clear Indx Indx
    
end


 % Failsafe: garder indices continus !!!
    u=unique(G);
    uini=unique(G);
    du=diff(u);
    x=find(du ~= 1);
    while ~isempty(x) % refait les indices de G
        d=find(G >= x(1)+2);
        G(d)=G(d)-1;
        u=unique(G);
        du=diff(u);
        x=find(du ~= 1);
    end
    Map=[unique(G), Map(:,2:3)];% refait indices dans Map

%%% Save the output; optional;
if SaveFile ==1
    save([FileBase '.CluRes.mat'], 'T', 'G', 'Map', 'Par');
end
