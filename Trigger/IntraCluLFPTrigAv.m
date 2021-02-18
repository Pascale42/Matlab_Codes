FileBase = '60920_4-7_2';

% Vm = MakeVm('6921_1-2','get');

LoadCluRes(FileBase);
load([FileBase '.CluRes.mat']);
T = round(T/16);
% [EvtRes, EvtClu, EvtLabels, Labels] = LoadEvt([FileBase '.evt.cel'],1250);

% Cell = EvtRes(EvtClu==find(strcmp(Labels,'cel')));
% CurrPlay = EvtRes(EvtClu~=find(strcmp(Labels,'cel')));
% 
% Cell = reshape(Cell,2,length(Cell)/2);
% 
% CurrPlay = reshape(CurrPlay,2,length(CurrPlay)/2)';
% 
% [gRes Ind] = SelectPeriods(Res,Cell,'d',1);
% gClu = Clu(Ind);
% 
% [gRes Ind] = SelectPeriods(gRes,CurrPlay,'d',0);
% gClu = gClu(Ind);


myclu=find(G==10);
Tclu=T(myclu);
Gclu=G(myclu);

Par=LoadPar([FileBase '.xml']);

eeg=LoadBinary([FileBase '.eeg'],20, Par.nChannels);
Eeg=eeg(min(Tclu):max(Tclu)); % for light purposes

TrigEeg = TriggeredAv(Eeg,30,30,Gclu,Tlclu); 


TrigEeg = squeeze(TrigEeg);
% 
% figure
% imagesc([-30:30]/1.25,[1:size(AvVm,2)],unity(AvVm)')
% for i=1:size(AvVm,2)
%     subplotfit(i,size(AvVm,2));
%     plot(AvVm(:

