%function out=CorrPower(filename,win)
% window in msec
function out=CorrPower(filename,varargin)

win = DefaultArgs(varargin, {20});
win = win *1250/5;
load([filename '.osc.mat']);
%pow(:,6) = 
% win = 100*1250/5;
[xc lag ] = XCovMatrix(pow(1:5:end,:), win);
out{1} = xc;
out{2} = lag;