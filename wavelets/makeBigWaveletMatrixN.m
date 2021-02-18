function [wave,period,scale,coi] = makeBigWaveletMatrixN(FileBase,varargin)

[Channels,s0,dj,normalize] = DefaultArgs(varargin,{[],-1,0.01,0});

Par = LoadPar([FileBase '.xml']);
sr = Par.lfpSampleRate;
dt = 1./sr;

if s0 == -1
    s0=2*dt;
end

Eeg = LoadBinary([FileBase '.eeg'],Channels(1));
pad=1;
nd = length(Eeg);
J1=ceil((log(nd*dt/s0)/log(2))/dj);

if isempty(Channels) 
    Channels = 1:Par.nChannels;
end


for i=1:length(Channels)
    disp(['Channel # ' num2str(Channels(i)) ' being processed, out of ' num2str(length(Channels))])
    Eeg = LoadBinary([FileBase '.eeg'], Channels(i));
    [wave,period,scale,coi] = b_wavelet_lin(Eeg,dt,pad,dj,s0,J1);
    save(['Wavelet' FileBase '_Ch' num2str(Channels(i)) '.mat'],'wave','period','scale','coi','-v7.3');
end

