%
% function OutArgs = CoherenceDensity(FileBase, State, fMode, Channels, ResCoef, FreqRange, Overwrite, SaveFig);
%
% State: 'THE' or 'SWS'
% SaveFig: 1 to save report
% fMode: 'compute' and 'display' modes
% Channels, ResCoef, FreqRange: Don't give if you want same channels as the
% ones computed in LoadFilEeg

function OutArgs = CoherenceDensity(FileBase, State, fMode, Channels, ResCoef, FreqRange, Overwrite, SaveFig);

switch fMode

    case 'compute'

        if exist([FileBase '.' mfilename State '.mat'])>0 & ~Overwrite
            sprintf('Already computed!');return
        end

        Par = LoadPar([FileBase '.xml']);

        % loading;
        STA = load([FileBase '.sts.' State]);
        

        % loads all channels and resample
        
        if ~isempty(Channels)
            EEG=[];
            for i=Channels
                %         for i=1:length(Par.AnatGrps(1,1).Channels)
                eeg = LoadBinary([FileBase '.eeg'], i, Par.nChannels);
                eeg = SelectPeriods(eeg(:),STA,'c',1);
                eeg = resample(eeg,1,ResCoef);
                EEG=[EEG eeg];
            end
            OutArgs.FreqRange = FreqRange;
            OutArgs.reSampleRate = Par.lfpSampleRate/ResCoef;
            OutArgs.Channels = Channels;
        else
            load([FileBase '.LoadFilEeg' State '.mat']);
            OutArgs.Channels = Data.Channels;
            EEG = Data.EEG;
            OutArgs.reSampleRate = Data.reSampleRate;
            ResCoef = Data.ResCoef;
            OutArgs.FreqRange = Data.FreqRange;
        end

        

        % for coherence & spectrogram
        if strcmp(State,'THE')==1; WinLengthSec=5; end
          if strcmp(State,'the')==1; WinLengthSec=5; end      

        if strcmp(State,'GAMTHE')==1 | strcmp(State,'GAMSWS')==1; WinLengthSec=0.5; end

        if strcmp(State,'SWS')==1; WinLengthSec=15; end
        
        if strcmp(State,'SO')==1; WinLengthSec=15; end
        
         if strcmp(State,'BIC')==1; WinLengthSec=10; end

        OutArgs.WinLengthSample = 2^round(log2(WinLengthSec*OutArgs.reSampleRate));
        OutArgs.nFFT = 2*OutArgs.WinLengthSample;

        % save parameters
        OutArgs.EEG = EEG;
        OutArgs.State = State;
        OutArgs.Par = Par;
        

        %         [y,f,phi]= mtchd(OutArgs.EEG,OutArgs.nFFT,OutArgs.reSampleRate,OutArgs.WinLengthSample,[],3,'linear',[],FreqRange);

        [y, f, phi, yerr, phierr, phloc, pow]=...
            mtptchd(OutArgs.EEG, [], [], OutArgs.nFFT,OutArgs.reSampleRate,OutArgs.WinLengthSample,[],3,'linear',[],OutArgs.FreqRange);
        % or use mtptchd without implementing the clusters # and time and
        % it retuns the yerr to plot the confidence interval
        % for this use confplot(f, y, yerr(:,:,:,1))

        OutArgs.mtchd.y = y;
        OutArgs.mtchd.f = f;
        OutArgs.mtchd.phi = phi;
        OutArgs.mtchd.yerr = yerr;
        OutArgs.mtchd.phierr = phierr;
        OutArgs.mtchd.phloc = phloc;
        OutArgs.mtchd.pow = pow;

        save([FileBase '.' mfilename State '.mat'], 'OutArgs');


    case 'display'
        load([FileBase '.' mfilename State '.mat']);

        nChannels = length(OutArgs.Channels);

        figure(24)
        for Ch1=1:nChannels
            for Ch2 = Ch1:nChannels
                subplot(nChannels, nChannels, Ch1 + (Ch2-1)*nChannels);
                if(Ch1==Ch2)
                    plot(OutArgs.mtchd.f,20*log10(abs(OutArgs.mtchd.y(:,Ch1,Ch2))+eps),'r'); grid on; 
                                        ylabel('psd (dB)'); xlabel('Frequency'); grid on; % power spectrum density
                else
                    confplot(OutArgs.mtchd.f,OutArgs.mtchd.y(:,Ch1,Ch2),OutArgs.mtchd.yerr(:,Ch1,Ch2,1)); grid on; 
%                     ylabel('chd (dB)'); grid on; % coherence density
                    subplot(nChannels, nChannels, Ch2 + (Ch1-1)*nChannels);
%                     plot(OutArgs.mtchd.f,unwrap(OutArgs.mtchd.phi(:,Ch1,Ch2)),'k'); grid on; 
                    confplot(OutArgs.mtchd.f,unwrap(OutArgs.mtchd.phi(:,Ch1,Ch2)),unwrap(OutArgs.mtchd.phierr(:,Ch1,Ch2)),unwrap(OutArgs.mtchd.phierr(:,Ch1,Ch2)),'k'); grid on;
                    % unwrap means that the phase will be continuously
                    % reprensented: it avoids jumps; do mod for cyclic phase
                    %                     ylabel('ph shift (rd)'); grid on; % phase shift
                end
            end
%             ForAllSubplots('set(gca, ''FontSize'', 4)')
        end

        if SaveFig > 0;
            reportfig(gcf,['CoherenceDensity_' OutArgs.State], 0, ['File ' OutArgs.Par.FileName ', State ' OutArgs.State ', ch ' num2str(OutArgs.Channels)],150);
        end
        
        case 'displayDeg'
        load([FileBase '.' mfilename State '.mat']);

        nChannels = length(OutArgs.Channels);

        figure(34)
        for Ch1=1:nChannels
            for Ch2 = Ch1:nChannels
                subplot(nChannels, nChannels, Ch1 + (Ch2-1)*nChannels);
                if(Ch1==Ch2)
                    plot(OutArgs.mtchd.f,20*log10(abs(OutArgs.mtchd.y(:,Ch1,Ch2))+eps),'r'); grid on; 
                                         ylabel('psd (dB)'); xlabel('Frequency'); grid on; axis tight% power spectrum density
                else
                    confplot(OutArgs.mtchd.f,OutArgs.mtchd.y(:,Ch1,Ch2),OutArgs.mtchd.yerr(:,Ch1,Ch2,1)); grid on; ylim([0 1])
                     ylabel('chd (dB)'); grid on; % coherence density
                    subplot(nChannels, nChannels, Ch2 + (Ch1-1)*nChannels);
%                     plot(OutArgs.mtchd.f,unwrap(OutArgs.mtchd.phi(:,Ch1,Ch2)),'k'); grid on; 
                    confplot(OutArgs.mtchd.f,unwrap(OutArgs.mtchd.phi(:,Ch1,Ch2))*180/pi,unwrap(OutArgs.mtchd.phierr(:,Ch1,Ch2))*180/pi,unwrap(OutArgs.mtchd.phierr(:,Ch1,Ch2))*180/pi,'k'); grid on;
                    % unwrap means that the phase will be continuously
                    % reprensented: it avoids jumps; do mod for cyclic phase
                                         ylabel('ph shift (deg)'); grid on; % phase shift
                end
            end
%             ForAllSubplots('set(gca, ''FontSize'', 8)')
        end

        if SaveFig > 0;
            reportfig(gcf,['CoherenceDensity_' OutArgs.State], 0, ['File ' OutArgs.Par.FileName ', State ' OutArgs.State ', ch ' num2str(OutArgs.Channels)],150);
        end
        
        case 'displaymatrix'
            
        load([FileBase '.' mfilename State '.mat']);

        nChannels = length(OutArgs.Channels);
        
        % mean phase shift
        OutArgs.mtchd.phim=zeros(nChannels,nChannels);
        % mean coherence
        OutArgs.mtchd.ym=zeros(nChannels,nChannels);
        
        if strcmp(State, 'THE')  
            intv=[20:33];
        end
               if strcmp(State, 'SWS')  
            intv=[22:27];
        end
        
        for Ch1=1:nChannels
            for Ch2 = Ch1:nChannels
                if(Ch1 == Ch2)
                    OutArgs.mtchd.ym(Ch1,Ch2)=1;
                    OutArgs.mtchd.phim(Ch1,Ch2)=0;
                else
                    OutArgs.mtchd.ym(Ch1,Ch2)=mean(OutArgs.mtchd.y(intv,Ch1,Ch2));
                   OutArgs.mtchd.phim(Ch1,Ch2)=mean(unwrap(OutArgs.mtchd.phi(intv,Ch1,Ch2)));
                end
            end
        end
   
 figure(442)
 imagesc(OutArgs.mtchd.ym'); colorbar;coloration('dark');
 figure(542)
        for Ch1=1:nChannels
            for Ch2 = Ch1:nChannels
                subplot(nChannels, nChannels, Ch1 + (Ch2-1)*nChannels);
        
        polar([0 OutArgs.mtchd.phim(Ch1,Ch2)],[0 0.8],'-r'); 
            end
        end
        
                   clear Ch1 Ch2 intv   nChannels
        
        ymdiff=OutArgs_The.mtchd.ym - OutArgs_Sws.mtchd.ym;
         figure(642)
imagesc(ymdiff'); ;coloration('CorrMatDiff');  

        
end
