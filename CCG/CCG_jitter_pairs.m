%   function CCGJitterPairs = CCG_jitter_pairs(FileBase, CluPre,CluPost,varargin)
%
%  [BinSize,HalfBins,jscale,njitter,alpha, extraname] = DefaultArgs(varargin,{20,50,5,1000,0.01,'']});
%
%  Gets the CCG_jitter of all the pairs CluPre-CluPost
%  and sorts the significant excitatory and inhibitory pairs
%
function CCGJitterPairs = CCG_jitter_pairs(FileBase, CluPre,CluPost,varargin)

[BinSize,HalfBins,jscale,njitter,alpha, extraname] = DefaultArgs(varargin,{32,21,5,1000,0.01,''});

    Par=LoadPar([FileBase '.xml']); 
    BinSize = floor(Par.SampleRate/1000); %1 ms


tic;

nCells = length(CluPre);

GSPExcG=[];
GSPInhG=[];
CCGR=[];
ExcPairs=[];
InhPairs=[];

for n=1:nCells
  n
  [ccgR,tR,GSPE,GSPI,ccgJ]=CCG_Jitter(FileBase,CluPre(n,:),CluPost(n,:),'c',0,BinSize,HalfBins,jscale,njitter,alpha);
  GSPExcG=[GSPExcG; GSPE];
  GSPInhG =[GSPInhG; GSPI];
  CCGR= [CCGR ccgR];  % concatenates on the 2nd dimension: 81 x (2*nCells) x 2
  pre=CluPre(n,:);
  post=CluPost(n,:);
  if ~exist([FileBase '.CCGJitter'],'dir')
    s = sprintf('mkdir %s.CCGJitter', FileBase);
    system(s)
  end
  cd([FileBase '.CCGJitter'])
  save([FileBase '.CCGJitter_' num2str(CluPre(n,1)) '.'  num2str(CluPre(n,2)) '-' num2str(CluPost(n,1)) '.' num2str(CluPost(n,2)) '.mat'], 'ccgR', 'tR', 'GSPE','GSPI','ccgJ','pre','post');
  cd ../

  MonoWindow1=find((tR>=1)&(tR<=7));     % Mono Window 1ms ~ 3ms
  MonoWindow2=find((tR>=-7)&(tR<=-1));

  %%%%%%%%%%% Excitatory
  if     any(GSPE(MonoWindow1)==1)
    ExcPairs=[ExcPairs; CluPre(n,:),CluPost(n,:)];
  elseif any(GSPE(MonoWindow2)==1)
    ExcPairs=[ExcPairs; CluPost(n,:),CluPre(n,:)];
  end
  %%%%%%%%%%% Inhibitory
  if     any(GSPI(MonoWindow1)==1)
    InhPairs=[InhPairs; CluPre(n,:),CluPost(n,:)];
  elseif any(GSPI(MonoWindow2)==1)
    InhPairs=[InhPairs; CluPost(n,:),CluPre(n,:)];
  end
end


%%%%%%%%%%%% Concatenate and save
CCGJitterPairs.ExcPairs=ExcPairs;
CCGJitterPairs.InhPairs=InhPairs;
CCGJitterPairs.GSPExcG=GSPExcG;
CCGJitterPairs.GSPInhG=GSPInhG;
CCGJitterPairs.ccgR=CCGR;
CCGJitterPairs.tR=tR;
CCGJitterPairs.CluPre=CluPre;
CCGJitterPairs.CluPost=CluPost;
CCGJitterPairs.njitter=njitter;
CCGJitterPairs.jscale=jscale;
CCGJitterPairs.alpha=alpha;

save([FileBase '.CCGJitterPairs' extraname '.mat'], 'CCGJitterPairs');
endt=toc

% 
%      load([FileBase '.CCGJitterPairs.mat']);
%      figure(45000)
% for j=1:size(CluPre,1)
%          bar(tR,ccgR(:,1,2),'w')
%       line(tR,ccgJ.m,'linestyle','--','color','b')
%       line(tR,ccgJ.ptMax,'linestyle','--','color','r')
%       line(tR,ccgJ.gbMax,'linestyle','--','color','m')
%       line(tR,ccgJ.ptMin,'linestyle','--','color','r')
%       line(tR,ccgJ.gbMin,'linestyle','--','color','m')
%       set(gca,'XLim',[min(tR),max(tR)])
%       title(['ccg   ' num2str(Pre(j,:))  ' - ' num2str(Post(j,:))]);
%       waitforbuttonpress
%     end
