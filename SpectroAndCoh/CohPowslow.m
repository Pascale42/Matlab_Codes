function out=CohPowslow(filename,varargin)
load([filename '.osc.mat']);
[y0 f0] =mtchd(pow(1:25:end,:),2^13,50);
out{1}=y0;
out{2}=f0;