%function PhaseDiscretize
filename = '836604ws';

FreqRange =[10 20];
EegCh=2;
WhatPhase =pi;
Thresh=3;
epsilon =0.05;
Par = LoadPar([filename '.par']);
SampleRate=round(1e6/Par.SampleTime);
eegSampleRate = 1250;

Eeg=readsinglech([filename '.eeg'],Par.nChannels,EegCh);

[WavePhase, WaveAmp] = WavePhase(Eeg, FreqRange);

GoodWavePeriods=find(((WaveAmp-mean(WaveAmp))/std(WaveAmp))>Thresh);

WhatPhaseT=find(abs(abs(WavePhase)-WhatPhase)<epsilon);

LockTimes=intersect(WhatPhaseT,GoodWavePeriods)*SampleRate/eegSampleRate;

sw =load([filename '.sw7']);
sw=sw(:,2);

t=[sw;LockTimes];
g=[ones(length(sw),1);2*ones(length(LockTimes),1)];
CCG(t,g,20*10,100,20000,[],'scale');