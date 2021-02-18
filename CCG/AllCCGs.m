%
% function AllCCGs(FileBase, where, mode, CluPre)
%
% where = 'ec' or 'hpc' or '8shk' or 'get'
% mode  = 'c' compute
%          d1 for figure display : recursive
%          d2 for display 1 clu vs all others;
%          d3 for d2 recursive
%


function [ccg, tbin, Map]=AllCCGs(FileBase, where, mode, CluPre);

Par = LoadPar([FileBase '.xml']);
 
 
if strcmp(where, 'ec') == 1
    spkgp=5:Par.nElecGps;
end


if strcmp(where, 'hpc') == 1
    spkgp=1:4;
end

if strcmp(where, '8shk') == 1
    spkgp=1:8;
end




if strcmp(mode,'c') == 1
    
    tic
    
    if strcmp(where, 'get') == 1
        spkgp= input('Enter SpikeGroups: ');
    end
    
    [T,G,Map]=LoadCluRes(FileBase,spkgp);
    
    BinSize = floor(Par.SampleRate/1000); % 1 ms
    HalfBins = 21; %  ? ? 20 ms
    
    [ccg, tbin] = CCG(T, G, BinSize, HalfBins, Par.SampleRate);
    
    save([FileBase '.AllCCGs.' where '.mat'], 'ccg', 'tbin', 'Map');
    toc
    
end



if strcmp(mode,'d1') == 1
    
    load([FileBase '.AllCCGs.' where '.mat']);
    
    x=ceil(sqrt(size(ccg,3)));
    
    figure(1)
    
    for n=1:size(ccg,3)
        for m=1:size(ccg,3)
            subplot(x,x,m)
            if n ==m 
            bar(tbin,ccg(:,n,m),'r')
            axis tight
            title([num2str(n) ' vs '  num2str(m)])
            else
                bar(tbin,ccg(:,n,m))
            axis tight
            title([num2str(n) ' vs '  num2str(m)])
            end
        end
        waitforbuttonpress
    end
    
    
end


if strcmp(mode,'d2') == 1
    
    load([FileBase '.AllCCGs.' where '.mat']);
    
    [~,~, idElCluPre] = Intersection(CluPre, Map(:,2:3)); clear dum;
    
    x=ceil(sqrt(size(ccg,3)));
    
    figure
    
    
    for m=1:size(ccg,3)
        subplot(x,x,m)
        bar(tbin,ccg(:,idElCluPre,m))
        axis tight
        title([num2str(idElCluPre) ' vs '  num2str(m)])
    end
    
    
end

if strcmp(mode,'d3') == 1  % recursive of 2
    
    load([FileBase '.AllCCGs.' where '.mat']);
    
    x=ceil(sqrt(size(ccg,3)));
    
    figure(1)
    
    for n=1:size(CluPre,1)
        clu=CluPre(n,:);
        [~,~, idElCluPre] = Intersection(clu, Map(:,2:3)); clear dum;
        for m=1:size(ccg,3)
            subplot(x,x,m)
            bar(tbin,ccg(:,idElCluPre,m))
            axis tight
            title([num2str(idElCluPre) ' vs '  num2str(m)])
        end
        waitforbuttonpress
    end
end


if strcmp(mode,'d4') == 1  % tout
    
    load([FileBase '.AllCCGs.' where '.mat']);
    
    figure
nclu=size(ccg,3);
%         for m=1:size(ccg,3)
%             subplot(23,23,m)
%             bar(tbin,ccg(:,idElCluPre,m))
%             axis tight
%             title([num2str(idElCluPre) ' vs '  num2str(m)])
%         end

        for clu1=1:nclu
            for clu2 = clu1:nclu
                subplot(nclu, nclu, clu1 + (clu2-1)*nclu);
                if(clu1==clu2)
                    bar(tbin,ccg(:,clu1,clu2),'r'); axis tight; title([num2str(clu1)]);set(gca, 'XTickLabel', {})
                else
                    bar(tbin,ccg(:,clu1,clu2)); axis tight;set(gca, 'XTickLabel', {})
                end
            end
            ForAllSubplots('set(gca, ''FontSize'', 6)')
        end

end





%%
%     Tclu=T(find(G==i));
%     Tclu = Tclu/20000; % time in sec
%     FreqClu=1*length(Tclu)/max(Tclu); % get frequency in Hz
%
%     legend(['Freq = ' num2str(FreqClu) ' Hz'])
%
%     waitforbuttonpress
%
%     ct=input('Cell type : Pyr=1  IN=2  Not Clear=3    ==> ')
%     if isempty(ct)
%         ct=input('Cell type: Pyr=1  IN=2  Not Clear=3    ==> ')
%     end
%     CT=[CT;ct];
%
% end



%% CCG of all the cells vs all the cells
% figure(1)
% for i=1:length(UG)
%     for j = 1:length(UG)
%     for n=1:length(j)
%         subplotfit(n,length(j))
%             %subplot(1,2,1)
%         bar(tbin,ccg(:,find(UG==j(n)),i))
%         axis tight
%         title(['cell pairs ' num2str(j(n)) ' vs ' num2str(UG(i))])
%     end
%
%     waitforbuttonpress
%     end
% end
