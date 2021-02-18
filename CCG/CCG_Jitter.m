%
%  [ccgR,tR,GSPExc,GSPInh,ccgJ]=CCG_Jitter(FileBase,CluPre,CluPost,fMode,| Display, BinSize,HalfBins,jscale,njitter,alpha)
%
%   [Display, BinSize,HalfBins,jscale,njitter,alph] = ...
%   DefaultArgs(varargin,{0,[],40,5,1000,0.01});
%
%  ***** Computes (and displays) CCG Jitter of one pair *****
%
%  fMode: 'c'  'd'  'dr'
%  ccgR is given in Hz, not in counts, and is computed for a 32,556 Hz samplig rate
%  
%  fMode: 'c' compute and 'd' display 
%         'dr' for recursive display of a list of pairs
%   
%  CluPre and CluPost   : [2 3]  (for shank 2 clu3)
%  Bizsize, HalfBins : --> like in 'CCG.m'
%  jscale            : jittering scale, unit is 'ms'
%  njitter           : # of jittering
%  alpha             : significant level
%   
%  ---------------------------------------------------
%  ccgR   : ccg of real data   <-- [ccg,t]=CCG(...);
%  tR       : t of real data
%  GSPExc : Global Significant Period of Mono Excitation.
%  GSPInh : Global Significant Period of Mono Inhibition.
%  ---------------------------------------------------
%   
%
%  Coded by  Shigeyoshi Fujisawa
%  based on Asohan Amarasimghan's resampling methods
%  Modified by Pascale Quilichini
%   

function [ccgR,tR,GSPExc,GSPInh,ccgJ]=CCG_Jitter(FileBase,CluPre,CluPost,fMode, varargin)

[Display,BinSize,HalfBins,jscale,njitter,alpha] = ...
  DefaultArgs(varargin,{0,[],21,5,1000,0.01});

%%% was:  [ccgR,tR,GSPExc,GSPInh,ccgjMtx]=CCG_jitter(spiket,spikeind,clu1,clu2,BinSize,HalfBins,jscale,njitter,alpha,PresenBoolen)

switch fMode

  case 'c'   % compute 

    Par=LoadPar([FileBase '.xml']); 
    BinSize = floor(Par.SampleRate/1000); % 1 ms
    
    %%%%%% Loads Res and get CCG
    if FileExists([FileBase '.CluRes.mat'])
      load([FileBase '.CluRes.mat']);
    else
      [T,G,Map,Par]=LoadCluRes(FileBase, 1);
    end
    [~,~, idElCluPre] = Intersection(CluPre, Map(:,2:3));
    [~,~, idElCluPost] = Intersection(CluPost, Map(:,2:3));
    Tpre = T(G==idElCluPre);
    Tpost = T(G==idElCluPost);
    clear idElCluPre idElCluPost 

    [ccgR, tR] = CCG([Tpre;Tpost],[ones(size(Tpre));2*ones(size(Tpost))], BinSize, HalfBins, Par.SampleRate,[1,2],'hz');


    %%%%%%  CCG for jittering data
    for i=1:njitter
      Tpost_jitter = Tpost + 2*(20*jscale)*rand(size(Tpost))-1*20*jscale;
      [ccg, tJ] = CCG([Tpre;Tpost_jitter],[ones(size(Tpre));2*ones(size(Tpost))], BinSize, HalfBins, Par.SampleRate,[1,2],'hz');
      ccgj(:,i)=ccg(:,1,2);
      ccgjmax(i)=max(ccgj(:,i));
      ccgjmin(i)=min(ccgj(:,i));
    end


    %%%%%%  Computes the pointwise line
    signifpoint = njitter*alpha;
    for i=1:length(tJ)
      sortjitterDescend  = sort(ccgj(i,:),'descend');
      sortjitterAscend   = sort(ccgj(i,:),'ascend');
      ccgjptMax(i) = sortjitterDescend(signifpoint);
      ccgjptMin(i) = sortjitterAscend(signifpoint);
    end


    %%%%%%  Compute the global line
    sortgbDescend   = sort(ccgjmax,'descend');
    sortgbAscend    = sort(ccgjmin,'ascend');
    ccgjgbMax  = sortgbDescend(signifpoint)*ones(size(tJ));
    ccgjgbMin  = sortgbAscend(signifpoint)*ones(size(tJ));

    ccgjm  = mean(ccgj,2);


    %%%%%%%%%%%%%% Significant Period
    findExc = find((ccgR(:,1,2)>=ccgjgbMax')&(ccgR(:,1,2)>0));
    findInh = find((ccgR(:,1,2)<=ccgjgbMin')&(ccgjgbMin'>0));

    GSPExc = zeros(size(tR));  % Global Significant Period of Mono Excitation
    GSPInh = zeros(size(tR));  % Global Significant Period of Mono Inhibition

    GSPExc(findExc) = 1;
    GSPInh(findInh) = 1;

    %     ccgjMtx=ccgj;


    %%%%%%%%%%%%%%% Presentation
     
    if Display==1
        figure(5000)
      %plot(tR,s2s.ccg(:,1,2),'color','k')
      bar(tR,ccgR(:,1,2),'k')
      line(tR,ccgjm,'linestyle','--','color','b','LineWidth',1.5)
      line(tR,ccgjptMax,'linestyle','--','color','r','LineWidth',1.5)
      line(tR,ccgjgbMax,'linestyle','--','color','m','LineWidth',1.5)
      line(tR,ccgjptMin,'linestyle','--','color','r','LineWidth',1.5)
      line(tR,ccgjgbMin,'linestyle','--','color','m','LineWidth',1.5)
      line(zeros(2),[0.1 ccgjgbMax(end)],'linestyle','--','color','r')
      set(gca,'XLim',[min(tR),max(tR)])
      title(['ccg   ' num2str(CluPre)  ' - ' num2str(CluPost)]);
    end
    %%%%%%%%%%%%%%% 
    
    ccgJ.m=ccgjm;
    ccgJ.ptMax=ccgjptMax;
    ccgJ.gbMax=ccgjgbMax;
    ccgJ.ptMin=ccgjptMin;
    ccgJ.gbMin=ccgjgbMin;
    % save([FileBase '.CCGJitter.mat']);

    if ~exist([FileBase '.CCGJitter'],'dir')
      s = sprintf('mkdir %s.CCGJitter', FileBase);
      system(s);
    end
    dn=[FileBase '.CCGJitter'];
    cd(dn)
    save([FileBase '.CCGJitter_' num2str(CluPre(1,1)) '.'  num2str(CluPre(1,2)) '-' num2str(CluPost(1,1)) '.' num2str(CluPost(1,2)) '.mat'],'ccgJ','ccgR','tR')
    cd ..

    
    
    
  case 'd'    % display 
    
    dn=[FileBase '.CCGJitter'];
    cd(dn)
    if FileExists([FileBase '.CCGJitter_' num2str(CluPre(1,1)) '.'  num2str(CluPre(1,2)) '-' num2str(CluPost(1,1)) '.' num2str(CluPost(1,2)) '.mat']);
        
     load([FileBase '.CCGJitter_' num2str(CluPre(1,1)) '.'  num2str(CluPre(1,2)) '-' num2str(CluPost(1,1)) '.' num2str(CluPost(1,2)) '.mat']);
    cd ..

    figure(4200)
    bar(tR,ccgR(:,1,2),'k')
    line(tR,ccgJ.m,'linestyle','--','color','b')
    line(tR,ccgJ.ptMax,'linestyle','--','color','r','LineWidth',1.5)
    line(tR,ccgJ.gbMax,'linestyle','--','color','m','LineWidth',1.5)
    line(tR,ccgJ.ptMin,'linestyle','--','color','r','LineWidth',1.5)
    line(tR,ccgJ.gbMin,'linestyle','--','color','m','LineWidth',1.5)
    line(zeros(2),[0.1 ccgJ.gbMax(end)],'linestyle','--','color','r')
    set(gca,'XLim',[min(tR),max(tR)])
    title(['ccg   ' num2str(CluPre)  ' - ' num2str(CluPost)]);
    else
        cd ..
        disp('To compute')
        return
    end

  case 'dr'  % display recursive

    Pre=CluPre; Post=CluPost;
    for j=1:length(Pre(:,1))
      dn=[FileBase '.CCGJitter'];
      cd(dn)
      load([FileBase '.CCGJitter_' num2str(Pre(j,1)) '.'  num2str(Pre(j,2)) '-' num2str(Post(j,1)) '.' num2str(Post(j,2)) '.mat']);
      cd ..
      figure(45000)
      bar(tR,ccgR(:,1,2),'k')
      line(tR,ccgJ.m,'linestyle','--','color','b','LineWidth',1.5)
      line(tR,ccgJ.ptMax,'linestyle','--','color','r','LineWidth',1.5)
      line(tR,ccgJ.gbMax,'linestyle','--','color','m','LineWidth',1.5)
      line(tR,ccgJ.ptMin,'linestyle','--','color','r','LineWidth',1.5)
      line(tR,ccgJ.gbMin,'linestyle','--','color','m','LineWidth',1.5)
      line(zeros(2),[0.1 ccgJ.gbMax(end)],'linestyle','--','color','r')
      set(gca,'XLim',[min(tR),max(tR)])
      title(['ccg   ' num2str(Pre(j,:))  ' - ' num2str(Post(j,:))]);
      waitforbuttonpress
    end

end




%  Example  : [ccg,tbin] = CCG_jitter(spiket,spikeind,clu1,clu2,32,60,5,500,0.01);
%             --- binsize is 1ms (32khz), jitter time scale is 5ms, 500 times jittering, p<0.01

