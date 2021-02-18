%  function [ExcPairs,InhPairs,GSPExcG,GSPInhG] = CCG_jitter_group(spiket,spikeind,Cells,BinSize,HalfBins,jscale,njitter,alpha)
%
%
function [ExcPairs,InhPairs,GSPExcG,GSPInhG] = CCG_jitter_initial(FileBase, CluPre,CluPost,varargin)
[BinSize,HalfBins,jscale,njitter,alpha] = DefaultArgs(varargin,{20,50,5,1000,0.01});
tic;

nCells = length(CluPre);

ii=0;jj=0;

GSPExcG=[];
GSPInhG=[];
CCGR=[];

for i=1:nCells
  i
%   for j=1:nCells
%     if i<j
      [ccgR,tR,GSPE,GSPI]=CCG_jitter(FileBase,CluPre(i,:),CluPost(i,:),BinSize,HalfBins,jscale,njitter,alpha,0);
GSPExcG=[GSPExcG; GSPE];
GSPInhG =[GSPInhG; GSPI];
CCGR= [CCGR ccgR];  % concatenates on the 2nd dimension: 81 x (2*nCells) x 2
      
%       GSPExcG(i,j,:) = GSPE;
%       GSPInhG(i,j,:) = GSPI;

       MonoWindow1=find((tR>=1)&(tR<=7));     % Mono Window 1ms ~ 3ms
       MonoWindow2=find((tR>=-7)&(tR<=-1));
     
      %%%%%%%%%%% Excitatory
      if     any(GSPE(MonoWindow1)==1)
        ii=ii+1;
        ExcPairs(ii,1)=Cells(i);
        ExcPairs(ii,2)=Cells(j);
      elseif any(GSPE(MonoWindow2)==1)
        ii=ii+1;
        ExcPairs(ii,1)=Cells(j);
        ExcPairs(ii,2)=Cells(i);
      end
      %%%%%%%%%%% Inhibitory
      if     any(GSPI(MonoWindow1)==1)
        jj=jj+1;
        InhPairs(jj,1)=Cells(i);
        InhPairs(jj,2)=Cells(j);
      elseif any(GSPI(MonoWindow2)==1)
        jj=jj+1;
        InhPairs(jj,1)=Cells(j);
        InhPairs(jj,2)=Cells(i);
      end

%     end
%   end
end

endt=toc