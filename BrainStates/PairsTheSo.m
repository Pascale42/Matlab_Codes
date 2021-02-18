%
% function PairsTheSo(FileBase)

function PairsTheSo(FileBase)

load([FileBase '.CCGJ.GoodPairs.mat']);

for n=1:size(CnxCt.pairs,1)
    
    
    cd([FileBase '.CCGJitterSt'])
    load([FileBase '.CCGJitterSt_' num2str(CnxCt.pairs(n,1)) '.'  num2str(CnxCt.pairs(n,2)) '-' num2str(CnxCt.pairs(n,4)) '.' num2str(CnxCt.pairs(n,5)) '-T.mat']);
    load([FileBase '.CCGJitterSt_' num2str(CnxCt.pairs(n,1)) '.'  num2str(CnxCt.pairs(n,2)) '-' num2str(CnxCt.pairs(n,4)) '.' num2str(CnxCt.pairs(n,5)) '-S.mat']);
    cd ../
    
    dn=[FileBase '.CCGJitterOne'];
    cd(dn)
    load([FileBase '.CCGJitterOne_' num2str(CnxCt.pairs(n,1)) '.'  num2str(CnxCt.pairs(n,2)) '-' num2str(CnxCt.pairs(n,4)) '.' num2str(CnxCt.pairs(n,5)) '.mat']);
    cd ..
    
        
    figure(42)
    %%% Theta 
    subplot(1,7,1); bar(tbinT,ccgT(:,1,1),'FaceColor',[0 0 0.9],'EdgeColor',[0 0 0.9])
    title([num2str(CnxCt.spk.the(n,2)) ' Hz (' num2str(CnxCt.spk.the(n,1)) ')'])
    subplot(1,7,2); bar(tbinT,ccgT(:,1,2),'FaceColor',[0.5 0.5 0.5],'EdgeColor',[0.5 0.5 0.5]);
    line(tbinT,ccgJT.m,'linestyle','-','color','b')
    line(tbinT,ccgJT.ptMax,'linestyle','-','color','r')
    line(tbinT,ccgJT.gbMax,'linestyle','--','color','m')
    line(tbinT,ccgJT.ptMin,'linestyle','-','color','r')
    line(tbinT,ccgJT.gbMin,'linestyle','--','color','m')
    set(gca,'XLim',[min(tbinT),max(tbinT)])
    title(['ccg Theta']);
    subplot(1,7,3); bar(tbinT,ccgT(:,2,2),'FaceColor',[0 0 0.9],'EdgeColor',[0 0 0.9])
    title([num2str(CnxCt.spk.the(n,4)) ' Hz (' num2str(CnxCt.spk.the(n,3)) ')'])
    %%% Slow Oscil 
    subplot(1,7,4); bar(tbinS,ccgS(:,1,1),'FaceColor', [0 0.6 0],'EdgeColor',[0 0.6 0]);
    title([num2str(CnxCt.spk.so(n,2)) ' Hz (' num2str(CnxCt.spk.so(n,1)) ')'])
    subplot(1,7,5); bar(tbinS,ccgS(:,1,2),'FaceColor',[0.5 0.5 0.5],'EdgeColor',[0.5 0.5 0.5]);
    line(tbinS,ccgJS.m,'linestyle','-','color','b')
    line(tbinS,ccgJS.ptMax,'linestyle','-','color','r')
    line(tbinS,ccgJS.gbMax,'linestyle','--','color','m')
    line(tbinS,ccgJS.ptMin,'linestyle','-','color','r')
    line(tbinS,ccgJS.gbMin,'linestyle','--','color','m')
    set(gca,'XLim',[min(tbinS),max(tbinS)]);
    title(['ccg SO']);
    subplot(1,7,6); bar(tbinS,ccgS(:,2,2),'FaceColor', [0 0.6 0],'EdgeColor',[0 0.6 0]);
    title([num2str(CnxCt.spk.so(n,4)) ' Hz (' num2str(CnxCt.spk.so(n,3)) ')'])
    %%% Whole
    subplot(1,7,7);
    bar(tR,ccgR(:,1,2),'k')
    line(tR,ccgJ.m,'linestyle','--','color','b')
    line(tR,ccgJ.ptMax,'linestyle','--','color','r','LineWidth',1.5)
    line(tR,ccgJ.gbMax,'linestyle','--','color','m','LineWidth',1.5)
    line(tR,ccgJ.ptMin,'linestyle','--','color','r','LineWidth',1.5)
    line(tR,ccgJ.gbMin,'linestyle','--','color','m','LineWidth',1.5)
    line(zeros(2),[0.1 ccgJ.gbMax(end)],'linestyle','--','color','r')
    set(gca,'XLim',[min(tR),max(tR)])
    title(['ccg   ' num2str(CnxCt.pairs(n,[1 2]))  ' - ' num2str(CnxCt.pairs(n,[4 5])) ' (' num2str(n) ')']);
    
    waitforbuttonpress;
end
