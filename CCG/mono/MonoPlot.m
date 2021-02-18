function MonoPlot(FileBase, ElClu1, ElClu2, State)

Els = unique([ElClu1(:,1); ElClu2(:,1)]);
load([FileBase '.thpar.mat']);

Par = LoadPar([FileBase '.xml']);
MinPeriod = 1; %seconds for theta periods selection
%[ThPhAll, ThAmpAll, ThFrAll] = ThetaParams(FileBase);

Period = load([FileBase '.sts.' State]);
Period = Period(find(diff(Period,1,2)>Par.lfpSampleRate*MinPeriod),:);

ThPh = SelectPeriods(ThPh, Period,'c',1);

[Res Clu Map] = LoadCluRes(FileBase,Els);

[rRes ind]  = SelectPeriods(round(Res*Par.lfpSampleRate/Par.SampleRate),Period,'d',1,1);
Clu = Clu(ind);
Res = Res(ind); 
n=size(ElClu1,1);

for ii=1:n
    refclu = find(Map(:,2)==ElClu1(ii,1) & Map(:,3)==ElClu1(ii,2));
    xclu = find(Map(:,2)==ElClu2(ii,1) & Map(:,3)==ElClu2(ii,2));
    
    figure(3348);clf
%     h = TrigRasters(Res(Clu==refclu), Par.lfpSampleRate*500/1000, Res(Clu==xclu),ones(sum(Clu==xclu),1), Par.lfpSampleRate, 1, ThPh(Res(Clu==refclu)));
%     set(h,'MarkerSize',6);
   [TrLag TrInd TrClu] = TrigRasters(Res(Clu==refclu), Par.SampleRate*500/1000, Res(Clu==xclu),ones(sum(Clu==xclu),1), Par.SampleRate, 1, ThPh(rRes(Clu==refclu)));
   subplot(211)
    hist2([TrLag TrInd],40,20);
    subplot(212)
    plot(TrLag,TrInd,'.'); xlim([-30 30]);

  % keyboard
   % waitforbuttonpress
end
%  [TrLag TrInd TrClu] = TrigRasters(Res(Clu==refclu), Par.SampleRate*500/1000, ...
%      Res(Clu==xclu),ones(sum(Clu==xclu),1), Par.SampleRate, 1, ThPh(rRes(Clu==refclu)));
%  
%  figure
%    subplot(211)
%     hist2([TrLag TrInd],40,20);
%     subplot(212)
%     plot(TrLag,TrInd,'.'); xlim([-30 30]);
%  
% 
%     
% 
%  [TrLag TrInd TrClu] = TrigRasters(shRes3, 1250*500/1000, ...
%      rRes(Clu==xclu),ones(sum(Clu==xclu),1), 1250, 1, ThPh(shRes3));
%  figure
%  hist2([TrLag TrInd],40,20);
%  
 