%function Peaks = EEGPeak(filename, Chan, FreqRange, Sign, Thr)
function Peaks = EEGPeak(filename,varargin)

[Chan, FreqRange, Sign, Thr] = DefaultArgs(varargin, {[1], [7 17], -1, 5});
% csdChans = [Chans(1)-2:Chans(end)+2];
% csd = MakeCSD(filename, csdChans, 2, 0);
% csd = sum(csd,2);
Par = LoadPar([filename '.par']);
eeg = readsinglech([filename '.eeg'],Par.nChannels,Chan);
eeg = FirFilter(eeg, 20, FreqRange/625, 'bandpass');
%use the sdame eeg variable to save memory1!!!!!!!!
%eeg = WaveletRec(eeg,FreqRange,1250,6,1);
%eeg=filtereeg(eeg,FreqRange,[],1250,20);
eeg = NormMatrix(eeg);
refr = 1250/(FreqRange(2));
Peaks = LocalMinimaK(-Sign*eeg, refr,  -Thr);

figure
hist(eeg,1000);
title(num2str(length(Peaks)));
if (nargout<1)
    msave([filename '.eegp' num2str(Chan) '-' num2str(round(mean(FreqRange))) ], Peaks);
end

Peaks = [Peaks eeg(Peaks)];