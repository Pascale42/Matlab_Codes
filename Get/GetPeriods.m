
% getting segments

[EvtRes, EvtClu, EvtLabels, Labels] = LoadEvt([FileBase '.evt.bon'],1250);
CelClu = find(strcmp(Labels,'ok'));
CelEvtInd = find(EvtClu==CelClu);
CelTime = EvtRes(CelEvtInd);
CelTime = sort(CelTime);
if rem(length(CelTime),2)~=0
    error('!! odd number of cel events!!');
end
nCel = length(CelTime)/2;
CelTime = reshape(CelTime,2,nCel)';
Cel = CelTime;
clear EvtRes EvtRes EvtLabels Labels CelClu CelEvtInd CelTime nCel



% load eeg


eeg = LoadBinary([FileBase '.eeg'], [1 2], Par.nChannels);
eeg = SelectPeriods(eeg(:),Cel,'c',1);



