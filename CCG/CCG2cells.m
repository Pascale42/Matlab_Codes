%
% Matlab function to show Cross-Corelograms between 2 cells
%
% function CCG2cells(FileBase,CluPre,CluPost,varargin)
% 
% [BinSize,HalfBins] = ...
%   DefaultArgs(varargin,{20,40});
%
% CluPre and CluPost are vectors such as [Shk Clu]


function [t, ccg]=CCG2cells(FileBase,CluPre,CluPost,ElGps, varargin)

[BinSize,HalfBins] = ...
  DefaultArgs(varargin,{20,40});


if FileExists([FileBase '.CluRes.mat'])
  load([FileBase '.CluRes.mat']);
else
  [T,G,Map,Par]=LoadCluRes(FileBase, 1);
   [T,G,Map]=LoadCluRes(FileBase, ElGps,0);
end

[~,~, idElCluPre] = Intersection(CluPre, Map(:,2:3));
[~,~, idElCluPost] = Intersection(CluPost, Map(:,2:3));

Tpre = T(find(G==idElCluPre));
Tpost = T(find(G==idElCluPost));

[ccg, t] = CCG([Tpre;Tpost],[ones(size(Tpre));2*ones(size(Tpost))], BinSize, HalfBins, 20000,[1,2],'count');
bar(t,ccg(:,1,2))


figure(7562956)
bar(t,ccg(:,1,2))


