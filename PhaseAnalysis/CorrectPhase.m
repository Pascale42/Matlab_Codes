function PhNew = CorrectPhase(Ph,varargin)
%function PhNew = CorrectPhase(Ph)
% correct continuous phase vector (ussualy form hilbert transform)
% for nonuniformity.
%does the ecdf transformatiton of the phase 
%so thatt PhNew = 2*pi*ecdf(Ph)-pi
% if Ph is uniform distributed then PhNew = Ph

[f,x,ind] = myecdf(Ph);

[dummy indr] = sort(ind);

PhNew = 2*pi*f(indr)-pi;


