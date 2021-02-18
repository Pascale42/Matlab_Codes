function out = MakeUniformDistr(in,varargin)
%function out = MakeUniformDistr(in,a,b)
%does the ecdf transformatiton of the in
%so thatt  = (b-a)2*pi*ecdf(in)+awhos
% if in is uniform distributed then out=in
[a,b] = DefaultArgs(varargin,{min(in), max(in)});
[f,x,ind] = myecdf(in);

[dummy indr] = sort(ind);

out = (b-a)*f(indr)+a;

