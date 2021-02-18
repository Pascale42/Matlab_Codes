
function All_CCGs(FileBase, SpkGrps)

Par = LoadPar([FileBase '.xml']);
BinSize = floor(Par.SampleRate/1000); % 1 ms
HalfBins = 21; %  ? ? 20 ms

if nargin < 2
    [T, G, Map]=LoadCluRes(FileBase);
else
    [T, G, Map]=LoadCluRes(FileBase, SpkGrps);
end
[T,it]=unique(T); G=G(it); clear it % prevoir si double detection de spike

% THE state
THE=load([FileBase '.sts.THE']);
THE=round(THE*(Par.SampleRate/Par.lfpSampleRate));
[Tt, Ind] = SelPerDiscr(T, THE,1,1);
Gt = G(Ind); clear Ind
[uClu, ~, NewClu] = unique(Gt);
x= ismember(Map(:,1), uClu);
mapt = Map(x,:); clear x
[ccgt, tbin] = CCG(Tt, Gt, BinSize, HalfBins, Par.SampleRate);

% SWS state
SWS=load([FileBase '.sts.SWS']);
SWS=round(SWS*(Par.SampleRate/Par.lfpSampleRate));
[Ts, Ind] = SelPerDiscr(T, SWS,1,1);
Gs = G(Ind); clear Ind
[uClu, ~, NewClu] = unique(Gs);
x= ismember(Map(:,1), uClu);
maps = Map(x,:); clear x
[ccgs, tbin] = CCG(Ts, Gs, BinSize, HalfBins, Par.SampleRate);

save([FileBase '.All_CCGs.mat'], 'ccgt', 'tbin', 'mapt','ccgs','maps');