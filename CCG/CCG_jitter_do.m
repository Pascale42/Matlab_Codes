%
%  function  [ccgR,tR,GSPExc,GSPInh,ccgJ]=CCG_jitter_do(FileBase,Tpre,Tpost,BinSize,HalfBins,Sampling, varargin)
%  [jscale,njitter,alpha,Display] = DefaultArgs(varargin,{5,1000,0.01,1});
% 
%  Revised June 2013

function [ccgR,tR,GSPExc,GSPInh,ccgJ]=CCG_jitter_do(FileBase,Tpre,Tpost,BinSize,HalfBins,Sampling, varargin)

[jscale,njitter,alpha,Display] = ...
  DefaultArgs(varargin,{5,1000,0.01,1});



%% Get CCG


[ccgR, tR] = CCG([Tpre;Tpost],[ones(size(Tpre));2*ones(size(Tpost))], BinSize, HalfBins, Sampling,[1,2],'hz');


%%  CCG for jittering data

for i=1:njitter
  Tpost_jitter = Tpost + 2*(20*jscale)*rand(size(Tpost))-1*20*jscale;
  [ccg, tJ] = CCG([Tpre;Tpost_jitter],[ones(size(Tpre));2*ones(size(Tpost))], BinSize, HalfBins, Sampling,[1,2],'hz');
  ccgj(:,i)=ccg(:,1,2);
  ccgjmax(i)=max(ccgj(:,i));
  ccgjmin(i)=min(ccgj(:,i));
end


%% Computes the pointwise line

signifpoint = njitter*alpha;
for i=1:length(tJ)
  sortjitterDescend  = sort(ccgj(i,:),'descend');
  sortjitterAscend   = sort(ccgj(i,:),'ascend');
  ccgjptMax(i) = sortjitterDescend(signifpoint);
  ccgjptMin(i) = sortjitterAscend(signifpoint);
end


%%  Compute the global line

sortgbDescend   = sort(ccgjmax,'descend');
sortgbAscend    = sort(ccgjmin,'ascend');
ccgjgbMax  = sortgbDescend(signifpoint)*ones(size(tJ));
ccgjgbMin  = sortgbAscend(signifpoint)*ones(size(tJ));

ccgjm  = mean(ccgj,2);


%% Significant Period

findExc = find((ccgR(:,1,2)>=ccgjgbMax')&(ccgR(:,1,2)>0));
findInh = find((ccgR(:,1,2)<=ccgjgbMin')&(ccgjgbMin'>0));

GSPExc = zeros(size(tR));  % Global Significant Period of Mono Excitation
GSPInh = zeros(size(tR));  % Global Significant Period of Mono Inhibition

GSPExc(findExc) = 1;
GSPInh(findInh) = 1;

%     ccgjMtx=ccgj;
ccgJ.m=ccgjm;
ccgJ.ptMax=ccgjptMax;
ccgJ.gbMax=ccgjgbMax;
ccgJ.ptMin=ccgjptMin;
ccgJ.gbMin=ccgjgbMin;


end
