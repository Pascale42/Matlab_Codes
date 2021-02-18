% Function out = myCSD_old(FileBase, State, fMode, ChannelRef, FreqRange,Overwrite, SaveFig);
%
% Loads LoadFilEeg and Computes CSD
%
% State: THE, SWS, GAMTHE or SWSTHE;
% fMode; use 'compute' otherwise 'display' to plot CSD;
% Winfeeg; window to plot feeg, in sec; 200 for theta; 1000 for SWS; 150 for gamma
%
% WinCSD; window for CSD: 500 for theta; 1500 for SWS; 100 for gamma
% Out args are saved in a .m file FileBase.myCSDTHE.mat (or
% FileBase.myCSDTHEGamma.mat)

function OutCsd = myCSD_old(FileBase, State, fMode, ChannelRef, FreqRange,Overwrite, SaveFig);

switch fMode
    case 'compute'

        if exist([FileBase '.' mfilename State num2str(ChannelRef) '.mat'])>0 & ~Overwrite; return
        end

       % Load worspace from LoadFilEeg;
%         load([FileBase '.LoadFilEeg' State '.mat'], 'Data');
        
        Data.STA = load([FileBase '.sts.' State]);
        Data.Par = LoadPar([FileBase '.xml']);
        [eeg Data.OrInd] = LoadBinary([FileBase '.eeg'], ChannelRef, Data.Par.nChannels,Data.STA);
        Data.eeg = resample(eeg',1,10);
        Data.reSampleRate = Data.Par.lfpSampleRate/10;
        Data.feeg = ButFilter(Data.eeg, 2, FreqRange/(Data.reSampleRate/2),'bandpass');
        

        % for coherence & spectrogram
        if strcmp(State,'THE')==1; OutCsd.Winfeeg=200; OutCsd.WinCSD=500; end

        if strcmp(State,'GAMTHE')==1 | strcmp(State,'GAMSWS')==1; OutCsd.Winfeeg=100; OutCsd.WinCSD=100;end

        if strcmp(State,'SWS')==1; OutCsd.Winfeeg=1000; OutCsd.WinCSD=1500;end

        OutCsd.ReferenceChannel = ChannelRef;
        ChRefId = ChannelRef;
        % CSD

        NotCloserThan = 0.5*Data.reSampleRate;
        allm = LocalMinima(-Data.feeg,NotCloserThan,0); 
        figure
        hist(log(Data.feeg(allm)),100)
        OutCsd.glm = allm(log(Data.feeg(allm))>7.5 & log(Data.feeg(allm))<8.5);
        
%         for visualizing:
%         segs = GetSegs(eeg,glm-70,141,[]);
%         imagesc(segs')
%         plot(mean(segs,2))
%         figure(3)



%         plot(Data.feeg(1:OutCsd.Winfeeg,ChRefId)) % 1000 for sws
%         legend(num2str(OutCsd.ReferenceChannel))
%         %             legend(num2str((Data.Channels)'))
%         title(['Phase, (feeg) -  Select threshold (Channel ' num2str(OutCsd.ReferenceChannel) ') to detect ' num2str(State)]);
%         xlabel('Time');
%         ylabel('Voltage, relative');
%         [dummy OutCsd.lmthr] = ginput(1);

%         NotCloserThan = 0.5*Data.reSampleRate;
%         allm = LocalMinima(-Data.feeg,NotCloserThan,0);
%         figure
%         hist(Data.feeg(allm),100)
        
        %[dummy OutCsd.lmthr] = ginput(1);

%         OutCsd.lm = LocalMinima(Data.feeg(:,ChRefId),10,OutCsd.lmthr);
        [OutCsd.EegSegAv OutCsd.EegSegStd OutCsd.EegSTrange] = ...
            TriggeredAvM(FileBase,Data.OrInd(OutCsd.glm),OutCsd.WinCSD,Data.Par.lfpSampleRate, Data.Par.nChannels,2);

        save([FileBase '.' mfilename State num2str(OutCsd.ReferenceChannel) '.mat'], 'OutCsd');


    case 'display'

        load([FileBase '.LoadFilEeg' State '.mat'], 'Data');
        load([FileBase '.myCSD' State  num2str(ChannelRef) '.mat'], 'OutCsd');

        figure
        ColorRange=PlotCSD(OutCsd.EegSegAv(:,((Data.Par.AnatGrps(1).Channels)+1)),OutCsd.EegSTrange,[],2, [], 2)

        if SaveFig == 1
            reportfig(gcf,['CSD_' State], 0, ['File ' FileBase ', State ' State, ', Ref Ch Id ' num2str(OutCsd.ReferenceChannel)],150);
        end
end
