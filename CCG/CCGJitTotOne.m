% 
%  CCGJitTotOne(FileBase, Pre, Post)

function CCGJitTotOne(FileBase, Pre, Post)

%% Loading

load([FileBase '.CluRes.mat']);
Pre=[2 2];    %%%%%%%%
Post=[2 8];   %%%%%%%%
[dum,dum, idElCluPre] = Intersection(Pre, Map(:,2:3));
[dum,dum, idElCluPost] = Intersection(Post, Map(:,2:3));
clear dum;
Tpre = T(find(G==idElCluPre));
Tpost = T(find(G==idElCluPost));


%% ACG, CCG and Jitter for significance

BinSize=20; % ie 1 ms
HalfBins=500;

[ccg,tbin,GSPE,GSPI,ccgJ] = CCG_jitter_do(FileBase,Tpre,Tpost,HalfBins);

if ~exist([FileBase '.CCGJitter'],'dir')
  s = sprintf('mkdir %s.CCGJitter', FileBase);
  system(s)
end
cd([FileBase '.CCGJitter'])
Clu.Pre=Pre; Clu.Post=Post;
save([FileBase '.CCGJitter_' num2str(Pre(1,1)) '.'  num2str(Pre(1,2)) '-' num2str(Post(1,1)) '.' num2str(Post(1,2)) '.mat'], 'ccg', 'tbin', 'GSPE','GSPI','ccgJ','Clu');
cd ../


%% Plotting now

figure(1)
subplot(231); bar(tbin,ccg(:,1,1),'b')
subplot(232); bar(tbin,ccg(:,1,2)); colormap([0.5 0.5 0.5]);
line(tbin,ccgJ.m,'linestyle','--','color','b')
line(tbin,ccgJ.ptMax,'linestyle','--','color','r')
line(tbin,ccgJ.gbMax,'linestyle','--','color','m')
line(tbin,ccgJ.ptMin,'linestyle','--','color','r')
line(tbin,ccgJ.gbMin,'linestyle','--','color','m')
set(gca,'XLim',[min(tbin),max(tbin)])
subplot(233); bar(tbin,ccg(:,2,2),'b')
t=[find(tbin==-30):find(tbin==30)];
subplot(234); bar(tbin(t),ccg(t,1,1),'b')
subplot(235); bar(tbin(t),ccg(t,1,2),'w')
line(tbin(t),ccgJ.ptMax(t),'linestyle','--','color','r')
line(tbin(t),ccgJ.gbMax(t),'linestyle','--','color','m')
line(tbin(t),ccgJ.ptMin(t),'linestyle','--','color','r')
line(tbin(t),ccgJ.gbMin(t),'linestyle','--','color','m')
subplot(236); bar(tbin(t),ccg(t,2,2),'b')
ForAllSubplots('axis tight')
