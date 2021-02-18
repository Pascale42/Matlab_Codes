%function y =FirFilter(x, n, wn, flag)
% x - signal, n - order, wn - center of filter, flag - band
function y =FirFilter(x, varargin)
[n, wn, flag] = DefaultArgs(varargin, {20, 50/625,'low'});
FiltVect = fir1(n,wn, flag);
nCh = size(x,2);
%  if size(x,1) == nCh
%      x=x';
%  end

for i=1:nCh
    y(:,i)=Filter0(FiltVect, x(:,i));
end

