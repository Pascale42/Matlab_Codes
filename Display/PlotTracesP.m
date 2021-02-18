function PlotTracesP(FileBase,Channels,FreqRange,t1,t2)

Par = LoadPar([FileBase '.xml']);

eeg  = LoadBinary([FileBase '.eeg'], Channels, Par.nChannels);
eeg=eeg';
vm = MakeVm(FileBase,'get');
feeg = ButFilter(eeg(:,1), 2, FreqRange/(Par.lfpSampleRate/2),'bandpass');
fvm = ButFilter(vm, 2, FreqRange/(Par.lfpSampleRate/2),'bandpass');

subplot(511)
plot (feeg((Par.lfpSampleRate*t1:Par.lfpSampleRate*t2)))
subplot(512)
plot (eeg((Par.lfpSampleRate*t1:Par.lfpSampleRate*t2),1))
subplot(513)
plot (eeg((Par.lfpSampleRate*t1:Par.lfpSampleRate*t2),2))
subplot(514)
plot (vm((Par.lfpSampleRate*t1:Par.lfpSampleRate*t2),:))
subplot(515)
plot (fvm((Par.lfpSampleRate*t1:Par.lfpSampleRate*t2)))