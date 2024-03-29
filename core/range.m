function y = range(x)
%RANGE  The range is the difference between the maximum and minimum values. 
%   Y = RANGE(X) calculates the range of the input.
%   For matrices RANGE(X) is a vector containing the range for each column.


y = max(x) - min(x);
