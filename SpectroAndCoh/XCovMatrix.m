function [xc, lag] = XCovMatrix(x, varargin)

[maxlag, flag] = DefaultArgs(varargin, {100, 'coeff'});
n =size(x,2);
xc =[];
for i=1:n
    for j=1:n
      %  if i~=j
            [xc(:,i,j), lag] = Sxcov(x(:,i),x(:,j),maxlag, flag);
       % end
    end
end




