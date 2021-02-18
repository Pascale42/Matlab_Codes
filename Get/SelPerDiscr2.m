%
% function [T, indT, indP] = SelPerDiscr2(t, Periods)
%
% Inputs ::
%           t :: vector of event timestamps
%           Periods :: n x 2 matrix of begin/end epochs
%
% Outputs ::
%           T : times of events t within Periods
%           indT : indices of these T
%           indP : indices of these Periods

function [T, indT, indP]=SelPerDiscr2(t,Periods)

indT=[];
indP=[];

for n=1:size(Periods,1)
    id=find(t >= Periods(n,1) & t <= Periods(n,2));
    indT=[indT; id];
    indP=[indP; repmat(n, size(id,1),1)];
    clear id
end
T=t(indT);


