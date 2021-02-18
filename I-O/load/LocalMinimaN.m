%function Mins = LocalMinimaN(p,AmpThr,NotCloserThen, DerStep, MaxDerPow)
% gives the coordinates of local minima of surface given by matrix p
% 
function [Mins MinsVal] = LocalMinimaN(p,AmpThr,NotCloserThan, varargin)

nDim = ndims(p);
sz = size(p);
[DerStep,MaxDerPow] = DefaultArgs(varargin,{1,1});
DerStep = 2*DerStep+1;
if length(NotCloserThan)==1 NotCloserThan = ones(nDim,1)*NotCloserThan;end

BelowThreshLin = find(p<AmpThr & ~isnan(p));
BelowThreshCol = Ind2Sub(sz,BelowThreshLin); % coordinates of points we will work with in each dimension
nBelow = length(BelowThreshLin);

%now make all possible shifts to nearest neighbor points. for robustness
%need to do 2 point shifts
HyperCube = reshape(ones(DerStep^nDim,1),DerStep*ones(nDim,1)'); 
Shifts = Ind2Sub(size(HyperCube),find(HyperCube)) - DerStep+1;
ShiftsNum = sum(Shifts~=0,2);
[ShiftsNum si] = sort(ShiftsNum);
Shifts = Shifts(si,:);

% use only the shifts of degree below or equal to MaxDerPow
% degree = number of dimensions shifted 
Shifts = Shifts(ShiftsNum<=MaxDerPow& ShiftsNum>0,:);
ShiftsNum = ShiftsNum(ShiftsNum<=MaxDerPow & ShiftsNum>0);
nShifts = size(Shifts,1);

IsMinima = logical(ones(nBelow,1));
ValueUnshifted = p(BelowThreshLin);
for s=1:nShifts
   % now we compute the sign of partial derivatives (value change for each
   % shift)
   ShiftsPos = BelowThreshCol + repmat(Shifts(s,:),nBelow,1);
 
   %get the position of islands shores 
   % in = ismember(ShiftedIndCol, BelowThreshCol,'rows');
   % ShoresPos = ShiftesPos(in,:);
       
   %consider shifts inside the boundaries - we don't count borders as minima yet
   GoodShifts = find(~any(ShiftsPos > repmat(sz,size(ShiftsPos,1),1) | ShiftsPos<1,2));
   
   %also will not consider the shifts inot NaNs
   GoodShifts = GoodShifts(~isnan( p(Sub2Ind(sz, ShiftsPos(GoodShifts,:))) ));
   
   %pOnTheBorderIndLin = find(~in); 
   %gives indices of elements on the border of islands BelowThreshold for this shift
      
   ValueShifted = p(Sub2Ind(sz, ShiftsPos(GoodShifts,:))); 
   IsMinima(GoodShifts) = IsMinima(GoodShifts) & (ValueShifted > ValueUnshifted(GoodShifts));
    
end
Mins = BelowThreshCol(IsMinima, :);
MinsVal = p(BelowThreshLin(IsMinima));

%now remove the small close minima
nMins= size(Mins,1);
TooClose = logical(ones(nMins,nMins));
for d=1:nDim
    TooClose = TooClose & abs(repmat(Mins(:,d),1,nMins)-repmat(Mins(:,d)',nMins,1))<NotCloserThan(d);
end
%    TooClose = tril(TooClose,-1);
TooClose =TooClose -diag(diag(TooClose));

[TooClose1 TooClose2] = find(TooClose);
nTooClose=length(TooClose1);
IsLarger = sparse(TooClose1, TooClose2,MinsVal(TooClose1)>MinsVal(TooClose2),nMins,nMins);
Delete = unique(find(full(sum(IsLarger,2))>0));
Mins(Delete,:) = [];
MinsVal(Delete) = [];


return

     
% nit=1;
% [dummy SortDimSize] = sort(sz,'descend');
% while 1
    %cheap approximate way - wrong one too :)
%     nMins= size(Mins,1);
%     TooClose = logical(ones(nMins-1,1));
%    
%     %[SortedMins SortMinsOrder] = sort(Mins(:,SortDimSize(1)),'ascend');
%     [Mins SortMinsOrder] = sortrows(Mins,SortDimSize);
%     MinsVal = MinsVal(SortMinsOrder);
%     for d=SortDimSize
%         
%         TooClose = TooClose & (diff(Mins(:,d)) < NotCloserThan(d));
% %         myTooClose = unique([myTooClose(:); myTooClose(:)+1]);
% %         TooClose = TooClose(SortMinsOrder(myTooClose));
%     end
%     TooClose = find(TooClose);
%     fprintf('Iteration %d, %d close pairs\n',nit,length(TooClose));
%     nit=nit+1;
%     if isempty(TooClose)
%         break;
%     end
%     Vals = [MinsVal(TooClose) , MinsVal(TooClose+1)];
%     [dummy Offset] = max(Vals,[],2);
%     Delete = unique(TooClose + Offset -1);
%     Mins(Delete,:) = [];
%     MinsVal(Delete) = [];
%end



