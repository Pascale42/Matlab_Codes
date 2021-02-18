%function out = PopulationBursts(Res,Clu,DownDur,UpDur,SampleRate)
% returns structure with vectors BurstTime, FirstBinNum, DownLen, BurstDur,
% SpkInBurst, FirstBinRate
function out = PopulationBursts(Res,Clu,DownDur,UpDur,SampleRate)

%sort the Res
[Res ind] = sort(Res);
Clu = Clu(ind);

[Burst, BurstLen, SpkPos, OutOf, FromBurst] = SplitIntoBursts(Res, DownDur);

nBurst = length(Burst);
BurstBeg = Res(Burst);
BurstEnd = Res(Burst+BurstLen-1);
BurstDur = BurstEnd-BurstBeg;

SpkTimeInBurst = Res - Res(Burst(FromBurst));
BurstBool = zeros(length(Res),1);
BurstBool(Burst)=1;

[dummy LagInd ] = histcI(SpkTimeInBurst,linspace(0,UpDur,10));

Ind = [LagInd Clu FromBurst];
SelInd = ~BurstBool & LagInd>0;% & SpkTimeInBurst>;
Ind = Ind(SelInd,:);
cnt = Accumulate(Ind,1,[max(Ind(:,1:2)) nBurst] );

HighFire = sum(cnt>0,3)./nBurst;
out.FirstBinNum = sq(sum(cnt(1,:,:)>0,2));

out.DownLen = [0; Res(Burst(2:end))-Res(Burst(2:end)-1)];

out.FirstBinRate = sq(sum(cnt(1,:,:),2));
out.BurstTime = Res(Burst);
out.BurstDur = BurstDur;
out.SpkInBurst = BurstLen;
