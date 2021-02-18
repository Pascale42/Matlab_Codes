%   function J = CCG_jitter_pairs_P(Pairs, T, G, Map, BinSize, HalfBins, Par,state, FileBase)
%
%  Gets the CCG_jitter of all the pairs CluPre-CluPost
%
function J = CCG_jitter_pairs_P(Pairs, T, G, Map, BinSize, HalfBins, Par,state, FileBase)

njitter=1000;
jscale=5;
alpha=0.01;
signifpoint = njitter*alpha;

  % Prep
  J.cc=[];
  J.m=[];
  J.ptMax=[];
  J.ptMin=[];
  J.gbMax=[];
  J.gbMin=[];
  
  for p=1:size(Pairs, 1)
      disp(num2str(Pairs(p,:)))
      [~,~, idElCluPre] = Intersection(Pairs(p,1:2), Map(:,2:3));
      [~,~, idElCluPost] = Intersection(Pairs(p,3:4), Map(:,2:3));
      Tpre = T(G==idElCluPre);
      Tpost = T(G==idElCluPost);
      clear idElCluPre idElCluPost
      
      [cc, J.tbin] = CCG([Tpre;Tpost],[ones(size(Tpre));2*ones(size(Tpost))], BinSize, HalfBins, Par.SampleRate,[1,2],'count');
      J.cc=[J.cc, squeeze(cc(:,1,2))];
      
      %%%%%%  CCG for jittering data
      for i=1:njitter
          Tpost_jitter = Tpost + 2*(20*jscale)*rand(size(Tpost))-1*20*jscale;
          [cj, tJ] = CCG([Tpre;Tpost_jitter],[ones(size(Tpre));2*ones(size(Tpost))], BinSize, HalfBins, Par.SampleRate,[1,2],'hz');
          Cj(:,i)=squeeze(cj(:,1,2));
          Cjmax(i)=max(Cj(:,i));
          Cjmin(i)=min(Cj(:,i));
      end
      
      
      %%%%%%  Computes the pointwise line
      for i=1:length(tJ)
          sortjitterDescend  = sort(Cj(i,:),'descend');
          sortjitterAscend   = sort(Cj(i,:),'ascend');
          ptMax(i) = sortjitterDescend(signifpoint);
          ptMin(i) = sortjitterAscend(signifpoint);
      end
      J.ptMax=[J.ptMax; ptMax];
      J.ptMin=[J.ptMin; ptMin];
      
      %%%%%%  Compute the global line
      sortgbDescend   = sort(Cjmax,'descend');
      sortgbAscend    = sort(Cjmin,'ascend');
      J.gbMax  = [ J.gbMax; sortgbDescend(signifpoint)*ones(size(tJ))];
      J.gbMin  = [J.gbMin; sortgbAscend(signifpoint)*ones(size(tJ))];
      
      J.m  = [J.m, mean(Cj,2)];
      
      clear sortjitterDescend sortjitterAscend sortgbDescend sortgbAscend  tJ...
          Cjmax Cjmin ptMax ptMin cc 
      
  end
save([FileBase '.Jitter.' state '.mat'], 'J');