%
% function MyColorMaps
%
% my upgraded JET colormap


function MyColorMaps

ici=pwd;

cd ~/Documents/MATLAB/
load('MyColormaps','mycmap')
set(figure(216793),'Colormap',mycmap)
cd(ici)
