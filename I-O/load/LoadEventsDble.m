function t= LoadEventsDble(FileBase, name, tag)
%
% function t= LoadEventsDble(FileBase, name, tag)
% 
% 


[EvtRes, EvtTags, EvtLabels, Labels] = LoadEvt([FileBase '.' name '.evt'],1250);


Tags = find(strcmp(Labels,tag));
EvtInd = find(EvtTags==Tags);
t = EvtRes(EvtInd);
t = sort(t);

if rem(length(t),2)~=0
    error('!! odd number of  events!!');
end
nTag = length(t)/2;
t = reshape(t,2,nTag)';





end

