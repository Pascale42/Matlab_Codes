% function myWavelets(data, FileBase, VARARGIN: period, FreqResol, Wcoefs, SR)
% 
% 
% [1 5] Hz ; 0.2,0.01,80
% [1 10]Hz ; 0.1,0.01,80
% [1 20]Hz ; 0.05,0.01,80
% 
% 
% period: period of the smallest wavelet (2 ms, for 500 Hz max) % once called s0
% 
% FreqResol: will determine the frequency resolution of wavelets % once called dj
% 
% Wcoefs: Number of wavelets coefficients  % once called mis
% with 81: get down 'till 100 Hz 
%(sinon il te descend jusqu'? une p?riode de 2s, soit 0.5 Hz, et surtout avec une
% r?solution fine comme ?a laisse tomber les d?g?ts en temps de computation). 
% 
% Donc, si tu veux par exemple aller que jusqu'? 400 Hz, tu joues sur s0 = 0.0025;
% Si tu veux descendre que jusqu'?  200 Hz, pour trouver le "mis" qu'il
% te faut: ouvre la variable "scale" qui a ?t? cr??e par le programme de
% wavelet. Elle contient les p?riodes des wavelets. La premiere est s0,
% puis tu rep?re l'indice de la p?riode max qui t'int?resse (en
% l'occurrence pour 200 Hz = 0.005 ms, donc mis=31);
% Si tu veux une meilleure r?solution en fr?quence tu diminues le "dj".

function myWavelets(data, FileBase, varargin)


[period, FreqResol, Wcoefs, SR] = ...
    DefaultArgs(varargin,{0.002, 0.0001, 81, 1250});


dt = 1./SR;
nd = length(data);

pad = 1; 


[wave,period,scale,coi] = b_wavelet_lin(data,dt,pad,FreqResol,period,-1,-1,-1,Wcoefs);

freq = 1./scale(1:Wcoefs);

power = abs(wave) .* abs(wave);

time = (1:nd) ./ SR;

Wout.wave=wave;
Wout.period=period;
Wout.scale=scale;
Wout.coi=coi;
Wout.FreqResol=FreqResol; 
Wout.Wcoefs=Wcoefs;

figure(1)

pcolor(time,freq,power)


shading 'interp'


% save([FileBase '.' mfilename '.mat'], 'time','freq','power','Wout');
