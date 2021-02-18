%
% Multitaper Time-Frequency Cross-Spectrogram
% Compute Phasogram and Coherence between Fields
%
% function OutArgs = CrossSpectrCoherence(FileBase, fMode, Channels, ResCoef, FreqRange, Overwrite, SaveFig);
%
% SaveFig: 1 to save report
% fMode: 'compute' and 'display' modes
% Channels, ResCoef, FreqRange: Don't give if you want same channels as the
% ones computed in LoadFilEeg

function OutArgs = CrossSpectrCoherence(FileBase, fMode, Channels, ResCoef, FreqRange, Overwrite, SaveFig);

switch fMode

    case 'compute'

        if exist([FileBase '.' mfilename '.mat'])>0 & ~Overwrite
            sprintf('Already computed!');return
        end

        Par = LoadPar([FileBase '.xml']);

       
        % loads all channels and resample
        
%         if ~isempty(Channels)
            EEG=[];
            for i=Channels
                %         for i=1:length(Par.AnatGrps(1,1).Channels)
                eeg = LoadBinary([FileBase '.eeg'], i, Par.nChannels);
                eeg = resample(eeg,1,ResCoef);
                EEG=[EEG eeg];
            end
            OutArgs.FreqRange = FreqRange;
            OutArgs.reSampleRate = Par.lfpSampleRate/ResCoef;
            OutArgs.Channels = Channels;



        % for coherence & spectrogram
        WinLengthSec=5;
        OutArgs.WinLengthSample = 2^round(log2(WinLengthSec*OutArgs.reSampleRate));
        OutArgs.nFFT = 2*OutArgs.WinLengthSample;

        % save parameters
        OutArgs.EEG = EEG;
        OutArgs.Par = Par;
  
        [y, f, t, phi] = mtchglong(OutArgs.EEG, OutArgs.nFFT,OutArgs.reSampleRate,OutArgs.WinLengthSample,[],3,'linear',[], FreqRange);
        OutArgs.mtchglong.y = y;
        OutArgs.mtchglong.f = f;
        OutArgs.mtchglong.t = t;
        OutArgs.mtchglong.phi = phi;

        save([FileBase '.' mfilename '.mat'], 'OutArgs');


    case 'display'
        load([FileBase '.' mfilename '.mat']);

        nChannels = length(OutArgs.Channels);

        figure(32)

        newplot; % take abs, and use image to display results
        for Ch1=1:nChannels, for Ch2 = Ch1:nChannels
                subplot(nChannels, nChannels, Ch1 + (Ch2-1)*nChannels);
                if Ch1==Ch2
                    if length(OutArgs.mtchglong.t)==1
                        imagesc([0 1/OutArgs.mtchglong.f(2)],OutArgs.mtchglong.f,20*log10(abs(OutArgs.mtchglong.y(:,:,Ch1,Ch2))+eps)');
                        axis xy; colormap(jet);
                    else
                        imagesc(OutArgs.mtchglong.t,OutArgs.mtchglong.f,20*log10(abs(OutArgs.mtchglong.y(:,:,Ch1,Ch2))+eps)');
                        axis xy; colormap(jet); axis off; %colorbar;
                    end
                    %                 title(['Power specgram ' num2str(Ch1)]);
                else
                    %imagesc the coherogram
                    imagesc(OutArgs.mtchglong.t,OutArgs.mtchglong.f,(abs(OutArgs.mtchglong.y(:,:,Ch1,Ch2)))');
                    axis xy; colormap(jet);axis off; %colorbar;
                    %                 title(['Coherogram ' num2str(Ch1)]);


                    %display phasogram
                    subplot(nChannels, nChannels, Ch2 + (Ch1-1)*nChannels);
                    imagesc(OutArgs.mtchglong.t,OutArgs.mtchglong.f,squeeze(OutArgs.mtchglong.phi(:,:,Ch1,Ch2))');
                    axis xy; colormap(jet);axis off %colorbar;
                    %                 title(['Phasogram ' num2str(Ch1)]);
                end
            end;
         ForAllSubplots('set(gca, ''FontSize'', 4)')
        end;
        xlabel('Time')
        ylabel('Frequency')

        if SaveFig > 0;
            reportfig(gcf,['CrossSpectrCoherence,_' OutArgs.Par.FileName '-' num2str(OutArgs.FreqRange) ' Hz'], 0, ['File ' OutArgs.Par.FileName ', ' num2str(OutArgs.FreqRange) ' Hz'],300);
        end
end
