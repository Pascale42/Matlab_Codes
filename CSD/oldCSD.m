% Function  = oldCSD(FileBase, fMode, State, RefChannels, FreqRange, Overwrite, SaveFig);
%
% Former myCSD
% Loads LoadFilEeg and Computes CSD on 4 shanks electrode
%
% State: THE, SWS, GAMTHE or SWSTHE;
% fMode; use 'compute' otherwise 'display' to plot CSD;
% Winfeeg; window to plot feeg, in sec; 200 for theta; 1000 for SWS; 150 for gamma
%
% WinCSD; window for CSD: 500 for theta; 1500 for SWS; 100 for gamma


function Csd = oldCSD(FileBase, fMode, State, Region, RefChannels, FreqRange,Overwrite, SaveFig);

switch fMode
    
    case 'compute'
        
        if exist([FileBase '.' mfilename State num2str(ChannelRef) '.mat'])>0 & ~Overwrite; return; end
        
        Par = LoadPar([FileBase '.xml']);
        STS = load([FileBase '.sts.' State]);
        
        
        if strcmp(Region,'hpc')==1;
            
            
            [eeg OrInd] = LoadBinary([FileBase '.eeg'], RefChannels, Par.nChannels,STS);
            eeg = resample(eeg',1,1);
            OrInd=resample(OrInd,1,1);
            reSampleRate = Par.lfpSampleRate/1;  
            feeg = ButFilter(eeg, 2, FreqRange/(reSampleRate/2),'bandpass');
            
            
            % for coherence & spectrogram
            if strcmp(State,'THE')==1; Winfeeg=200; WinCSD=500; end
            
            if strcmp(State,'GAMTHE')==1 | strcmp(State,'GAMSWS')==1; Winfeeg=100; WinCSD=100;end
            
            if strcmp(State,'SWS')==1; Winfeeg=1000; WinCSD=1500;end
            
            % CSD
            
            NotCloserThan = 0.5*reSampleRate;
            
            tic
            allm1 = LocalMinima(-feeg(:,1),NotCloserThan,0); % no loop since the size of the matrix changes
            glm1 = allm1(log(feeg(allm1,1))>7.5 & log(feeg(allm1,1))<8.5);
            [CSD(1).EegSegAv CSD(1).EegSegStd CSD(1).EegSTrange] = ...
                TriggeredAvMP(FileBase,OrInd(glm1),WinCSD,Par.lfpSampleRate, Par.nChannels,Par.SpkGrps(1).Channels+1);
            toc
            
            tic
            allm2 = LocalMinima(-feeg(:,2),NotCloserThan,0);
            glm2 = allm2(log(feeg(allm2,2))>7.5 & log(feeg(allm2,2))<8.5);
            [CSD(2).EegSegAv CSD(2).EegSegStd CSD(2).EegSTrange] = ...
                TriggeredAvMP(FileBase,OrInd(glm2),WinCSD,Par.lfpSampleRate, Par.nChannels,Par.SpkGrps(2).Channels+1);
            toc
            
            tic
            allm3 = LocalMinima(-feeg(:,3),NotCloserThan,0);
            glm3 = allm3(log(feeg(allm3,3))>7.5 & log(feeg(allm3,3))<8.5);
            [CSD(3).EegSegAv CSD(3).EegSegStd CSD(3).EegSTrange] = ...
                TriggeredAvMP(FileBase,OrInd(glm3),WinCSD,Par.lfpSampleRate, Par.nChannels,Par.SpkGrps(3).Channels+1);
            toc
            
            tic
            allm4 = LocalMinima(-feeg(:,4),NotCloserThan,0);
            glm4 = allm4(log(feeg(allm4,4))>7.5 & log(feeg(allm4,4))<8.5);
            [CSD(4).EegSegAv CSD(4).EegSegStd CSD(4).EegSTrange] = ...
                TriggeredAvMP(FileBase,OrInd(glm4),WinCSD,Par.lfpSampleRate, Par.nChannels,Par.SpkGrps(4).Channels+1);
            toc
            
            
            save([FileBase '.' mfilename '-' Region num2str(ChannelRef) '-' State '.mat'], 'CSD');
        end
        
    case 'display'
        
        %         load([FileBase '.LoadFilEeg' State '.mat'], 'Data');
        Par = LoadPar([FileBase '.xml']);
        load([FileBase '.' mfilename '-' State  '-' num2str(ChannelRef) '.mat'], 'Csd');
        
        figure
        subplot(141); PlotCSD(CSD(1).EegSegAv,CSD(1).EegSTrange,[],2, [], 2);
        subplot(142); PlotCSD(CSD(2).EegSegAv,CSD(2).EegSTrange,[],2, [], 2);
        subplot(143); PlotCSD(CSD(3).EegSegAv,CSD(3).EegSTrange,[],2, [], 2);
        subplot(144); PlotCSD(CSD(4).EegSegAv,CSD(4).EegSTrange,[],2, [], 2);
        
        %         if SaveFig == 1
        %             reportfig(gcf,['CSD_' State], 0, ['File ' FileBase ', State ' State, ', Ref Ch Id ' num2str(Csd.ReferenceChannel)],150);
        %         end
end
