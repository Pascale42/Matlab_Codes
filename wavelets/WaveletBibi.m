% function [power, freq, time, mpow] = WaveletBibi(data, saveFile, VARARGIN:: sr, s0, mis, dj)

% s0 = 0.02; -> (50 Hz) taille de la plus grande ondelette
% dj = 0.001; resolution en frequence
% mis = 1981;  81 -> 10Hz; 981 -> 1 Hz; 1981 -> 0.5 Hz taille de la plus petite ondelette

function [power, freq, time, mpow] = WaveletBibi(data, saveFile, varargin)

[sr, s0, mis, dj] = ...
    DefaultArgs(varargin,{1250, 0.02, 1981, 0.001});

[wave,~,scale] = b_wavelet_lin(data,(1./sr),1,dj,s0,-1,-1,-1,mis);
freq = 1./scale(1:mis);
power = abs(wave) .* abs(wave);
time = (1:length(data)) ./ sr;

mpow=mean(power,2);

if ~isempty(saveFile)
save([saveFile '.mat'], 'time', 'data', 'freq', 'power', 'mpow')
end

figure
subplot(2,3,[1 2])
plot(time, data)
subplot(2,3,[4 5])
pcolor(time,freq,power); shading 'interp'; colorbar
subplot(2,3,[3 6])
plot(freq, 20*log10(abs(mpow)+eps)); xlabel('Frequency, Hz'); ylabel('Power, dB')





%% FROM MARIE :: tu copie/paste ca:

% dt = 1./Par.lfpSampleRate;
% nd = length(eeg);
% s0 = 0.002; c'est la periode de la plus petite wavelette (s0 = 0.002 == 2 ms, pour faire du 500 Hz max)
% dj = 0.0001;  ca va determiner la resolution des wavelets en frequence
% pad = 1; ca on en parlera plus tard on s'en fout
% mis = 81; c'est le nombre de wavelets coefficients qui t'interessent. 
% avec 81 tu descends jusqu'a 100 Hz (sinon il te descend jusqu'a une periode de 2s, soit 0.5 Hz, 
% et surtout avec une resolution fine comme ca laisse tomber les degats en temps de computation). 
% pour savoir comment changer ce nombre a ton gre sans tatonner, cf en dessous.

% [wave,period,scale,coi] = b_wavelet_lin(eeg,dt,pad,dj,s0,-1,-1,-1,mis);
% freq = 1./scale(1:mis);
% power = abs(wave) .* abs(wave);
% time = (1:nd) ./ Par.lfpSampleRate;
% pcolor(time,freq,power)
% shading 'interp'
% freq(end)
% 
% Tadaaaa !!!
% 
% 
% Donc, si tu veux par exemple aller que jusqu'? 400 Hz, tu joues sur s0 = 0.0025;
% Si tu veux descendre que jusqu'a  200 Hz, pour trouver le "mis" qu'il te faut: 
% ouvre la variable "scale" qui a ete creee par le programme de wavelet. 
% Elle contient les periodes des wavelets. La premiere est s0,
% puis tu repere l'indice de la periode max qui t'interesse (en l'occurrence pour 200 Hz = 0.005 ms, donc mis=31);
% Si tu veux une meilleure resolution en frequence tu diminues le "dj".
