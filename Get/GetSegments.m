function Seg = GetSegments(FileBase, evt, label)

%% function Seg = GetSegments(FileName, evt, label)
%
% evt =  'evt.xxx' or 'xxx.evt'
% label = tag to extract from the evt file
% SR = sampling rate
% 
% Then use y=SelectPeriods(x, Seg, 'c', 1); to extract the eeg between the
% Seg timestamps


%% getting segments

[EvtT, EvtP, EvtLabels, Labels] = LoadEvt([FileBase '.' evt],1250);
Points = find(strcmp(Labels,label));
EvtInd = find(EvtP==Points);
Time = EvtT(EvtInd);
Time = sort(Time);

if rem(length(Time),2)~=0
    disp('!! odd number of events!!');
end

nEvt = length(Time)/2;
Seg = reshape(Time,2,nEvt)';


