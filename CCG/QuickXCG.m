
% function QuickXCG(FileBase,CluPre,CluPost);


function QuickXCG(FileBase,CluPre,CluPost);



[T,G,Map]=LoadCluRes(FileBase);

Par=LoadPar([FileBase '.xml']);

BinSize = floor(Par.SampleRate/1000); %1 ms
HalfBins = 21; %  ? ? 20 ms



if nargin <3  % ACG
    
    reply=1;
    while reply == 1
        
        [~,~, idElCluPre] = Intersection(CluPre, Map(:,2:3)); clear dum;
        Tpre = T(find(G==idElCluPre)); clear idElCluPre
        
        
        [ccg, tbin] = CCG(Tpre, ones(size(Tpre)), BinSize, HalfBins, Par.SampleRate,unique(ones(size(Tpre))),'hz');
        figure(42000)
        bar(tbin,ccg(:,1,1),'k')
        axis tight
        title(['ACG   ' num2str(CluPre)]);  ylabel('Hz')
        
        reply = input('Do you want more? 0/1: ');
        if reply == 1
            CluPre = input('Cluster: ');
        end
    end
    
else  % CCG
    
    reply=1;
    while reply == 1
        
        [~,~, idElCluPre] = Intersection(CluPre, Map(:,2:3));
        [~,~, idElCluPost] = Intersection(CluPost, Map(:,2:3)); clear dum
        Tpre = T(G==idElCluPre);
        Tpost = T(G==idElCluPost); clear idElCluPre idElCluPost
        
        
        [ccgR, tR] = CCG([Tpre;Tpost],[ones(size(Tpre));2*ones(size(Tpost))], BinSize, HalfBins, Par.SampleRate,[1,2],'hz');
        
        figure(42001)
        bar(tR,ccgR(:,1,2),'k')
        axis tight
        title(['CCG   ' num2str(CluPre)  ' - ' num2str(CluPost)]); ylabel('Hz')
        
        reply = input('Do you want more? 0/1: ');
        if reply == 1
            CluPre = input('Cluster Pre: ');
            CluPost= input('Cluster Post: ');
        end
    end
end








