cdk
FileBase='70720_1-2';
ChannelRef= 12;

cdp
FileBase='20110912PCA-2';
ChannelRef= 54;



cd(FileBase)
State='THE';
OutCsd.Winfeeg=200;
FreqRange=[2 6];
OutCsd.WinCSD=500;


Par = LoadPar([FileBase '.xml']);
STS = load([FileBase '.sts.' State]);

eeg = LoadBinary([FileBase '.eeg'], ChannelRef, Par.nChannels);
% eeg = resample(eeg,1,10);
[eeg orind]= SelectPeriods(eeg(:),STS,'c',1);

reSampleRate = Par.lfpSampleRate/10;
feeg = ButFilter(eeg, 2, FreqRange/(reSampleRate/2),'bandpass');


OutCsd.ReferenceChannel = ChannelRef;
ChRefId = ChannelRef;


NotCloserThan = 0.5*reSampleRate;
allm = LocalMinima(-feeg,NotCloserThan,0);
figure
hist(log(feeg(allm)),100)
OutCsd.glm = allm(log(feeg(allm))>6 & log(feeg(allm))<8);


%

%         for visualizing:
segs = GetSegs(eeg,OutCsd.glm-70,141,[]);
imagesc(segs'); hold on
plot(mean(segs,2))

figure
plot(feeg(1:OutCsd.Winfeeg)) 
legend(num2str(ChannelRef))
title(['Phase, (feeg) -  Select threshold (Channel ' num2str(ChannelRef) ') to detect ' num2str(State)]);
xlabel('Time');
ylabel('Voltage, relative');
[dummy OutCsd.lmthr] = ginput(1);

NotCloserThan = 0.5*reSampleRate;
allm = LocalMinima(-feeg,NotCloserThan,0);
figure
hist(feeg(allm),100)

[dummy OutCsd.lmthr] = ginput(1);

OutCsd.lm = LocalMinima(feeg(:),10,OutCsd.lmthr);
%
save('CSD-Test.mat')

[OutCsd.EegSegAv OutCsd.EegSegStd OutCsd.EegSTrange] = ...
    TriggeredAvM(FileBase,orind(OutCsd.glm),OutCsd.WinCSD,Par.lfpSampleRate,Par.nChannels,3, Par.AnatGrps(6).Channels+1);

[OutCsd.EegSegAv OutCsd.EegSegStd OutCsd.EegSTrange] = ...
    TriggeredAvM(FileBase,orind(OutCsd.glm),OutCsd.WinCSD,Par.lfpSampleRate,Par.nChannels,2);

figure
ColorRange=PlotCSD(OutCsd.EegSegAv(:,((Par.AnatGrps(6).Channels)+1)),OutCsd.EegSTrange,[],2, [], 2)