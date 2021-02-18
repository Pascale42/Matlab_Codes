% function Out=GetTags(FileName,lab,shape)
%
% shape is 1 or 2 columns

function Out=GetTags(FileName,lab, shape)


[EvtRes, Evtnum, EvtLabels, Labels] = LoadEvt(FileName,1250);
tag = find(strcmp(Labels,lab));
EvtInd = find(Evtnum==tag); Out = EvtRes(EvtInd); Out = sort(Out);

if shape==2
if rem(length(Out),2)~=0
    error('!! odd number of cel events!!'); else disp('----->  ok   :)')
end
nout = length(Out)/2; Out = reshape(Out,2,nout)'; 
end
clear EvtRes Evtnum EvtLabels Labels tag EvtInd nout




















% [EvtRes, Evtnum, EvtLabels, Labels] = LoadEvt('Francis-20100406-1to9.evt.rip',1250);
% tag = find(strcmp(Labels,'pic'));
% EvtInd = find(Evtnum==tag); Rip = EvtRes(EvtInd); Rip = sort(Rip);
% if rem(length(Rip),2)~=0
%     error('!! odd number of cel events!!'); else disp('----->  ok   :)')
% end
% nrip = length(Rip)/2; Rip = reshape(Rip,2,nrip)'; clear EvtRes Evtnum EvtLabels Labels tag EvtInd nrip