% Crossings = SchmittTrigger(Signal,UpThresh,DownThresh)
%
% finds all those points where the signal crosses UpThresh in the
% upwards direction, but only for the first time after each time it was below DownThresh

function [Crossings, UpCrossings, DownCrossings] = SchmittTrigger(Signal,UpThresh,DownThresh)

% this is the sort of thing that would be so much easier not vectorized

if numel(Signal)~=max(size(Signal)) % before :: prod(size(Signal))~=max(size(Signal))
	error('Can only take a vector input');
end
Signal = Signal(:);

n = length(Signal);
PrevVal = [(UpThresh+DownThresh)/2; Signal(1:n-1)];

UpCrossings = find(PrevVal<UpThresh & Signal>=UpThresh);
DownCrossings = find(PrevVal>DownThresh & Signal<=DownThresh);

% now sort the crossings in order of occurrence
UpDownUnsort = [UpCrossings; DownCrossings];
TypeUnsort = [ones(size(UpCrossings)) ; -ones(size(DownCrossings))];

% sort them
[UpDownSort, Index] = sort(UpDownUnsort);
TypeSort = TypeUnsort(Index);

% we want all those times a upcrossing comes after a downcrossing
PrevType = [NaN; TypeSort(1:length(TypeSort)-1)];

% Pascale 10/28/15


Crossings = UpDownSort(find(PrevType==-1 & TypeSort==1));


