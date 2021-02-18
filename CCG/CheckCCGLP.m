%
% function CheckCCGLP(FileBase,sg1, clu1, sg2, varargin)
%
% [BinSize, HalfBins] = DefaultArgs(varargin,{20,40});
%



function CheckCCGLP(FileBase,sg1, clu1, sg2, varargin)

[plotperc, plotccgs, BinSize, HalfBins] = DefaultArgs(varargin,{1,0,20,40});



if FileExists([FileBase '.CheckCCGLP-' num2str(sg1) '.' num2str(clu1) '-' num2str(sg2) '.mat'])
    load([FileBase '.CheckCCGLP-' num2str(sg1) '.' num2str(clu1) '-' num2str(sg2) '.mat']);
    disp('CCGs already computed')
else
    
    Par = LoadPar([FileBase '.xml']);
    
    
    
    % 2nd spikegroup to load
    g2=LoadClu([FileBase '.clu.' num2str(sg2)]);
    t2=load([FileBase '.res.' num2str(sg2)]);
    
    indx=find(g2>1); % remove clu 0 and 1
    g2=g2(indx);
    t2=t2(indx);
    %     clus2=length(unique(g2));
    
    % 1st spikegroup and clu1 to load
    g1=LoadClu([FileBase '.clu.' num2str(sg1)]);
    t1=load([FileBase '.res.' num2str(sg1)]);
    id1=find(g1 == clu1);
    t1=t1(id1);
    
    T=[t1; t2];
    G=[ones(size(t1)); g2];
    
    disp('computing CCGs: may take time...')
    [ccg, t] = CCG(T,G, BinSize, HalfBins, Par.SampleRate,unique(G),'count');
    
    save([FileBase '.CheckCCGLP-' num2str(sg1) '.' num2str(clu1) '-' num2str(sg2) '.mat'], 'ccg', 't')
    
end


if plotperc == 1
    
    
    %     Totccg=squeeze(sum(ccg(:,1,[2:end])));
    %     Pikccg=squeeze(ccg(HalfBins+1,1,[2:end]));
    
    Perccg=round((squeeze(ccg(HalfBins+1,1,[2:end]))).*100./(squeeze(sum(ccg(:,1,[2:end])))));
    
    MyColorMaps
    imagesc([2:size(ccg,3)],0.5,Perccg'); colorbar; caxis([0 100]);
    title('Click on the neurons for displaying CCGs')
    
    [x,dum]=ginput;
    x=round(x); clear dum
    
    if length(x) > 0
        
        for n=1:length(x)
            
            totccg=sum(ccg(:,1,x(n)));
            pikccg=ccg(HalfBins+1,1,x(n));
            perccg=round(pikccg*100/totccg);
            
            bar(t,ccg(:,1,x(n)));   axis tight
            title(['Central Bin is ' num2str(pikccg) '/' num2str(totccg) ' counts for ' num2str(sg1) '.' num2str(clu1) ' vs ' num2str(sg2) ...
                '.' num2str(x(n)) ' = ' num2str(perccg) '%'])
            waitforbuttonpress
            
        end
    end
    
    disp('Clusters to Check: ')
    x
end



% if plotccgs == 1
%
% figure(7562956);
%
% for n=1:clus2
%
% totccg=sum(ccg(:,1,n+1));
% pikccg=ccg(HalfBins+1,1,n+1);
% perccg=round(pikccg*100/totccg);
%
%
% bar(t,ccg(:,1,n+1));   axis tight
% title(['Central Bin is ' num2str(pikccg) '/' num2str(totccg) ' counts for ' num2str(sg1) '.' num2str(clu1) ' vs ' num2str(sg2) ...
%     '.' num2str(n+1) ' = ' num2str(perccg) '%'])
% waitforbuttonpress
%
%
% end
% end



% sg1=7; sg2=8; clu1=16; clu2=12;
% BinSize=20; HalfBins=40;
