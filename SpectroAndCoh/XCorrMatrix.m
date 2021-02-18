%function [xc, lag] = XCorrMatrix(x, varargin)
function [xc, lag] = XCorrMatrix(x, varargin)
[maxlag, flag] = DefaultArgs(varargin, {100, 'unbiased'});
n =size(x,2);
xc =[];
for i=1:n
    for j=1:n
        if i~=j
            [xc(:,i,j), lag] = xcorr(x(:,i),x(:,j),maxlag, flag);
        end
    end
end




