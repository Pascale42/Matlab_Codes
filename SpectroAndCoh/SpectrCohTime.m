% function [y, f, t, phi] = SpectrCohTime(eeg, fMode, VARARGIN: FreqRange, WinLengthSec, SR, Whiten, name)
%
% fMode:  'c'  :: compute and save
%              'd'  :: display; only if saved before
%             'cd' ::compute and display, not saved
%
% DefaultArgs: FreqRange=[0.1 50]; WinLengthSec=4;
%              SR(sampling rate)=1250; Whiten=1; name=[]; 
%
% Use 'name' to save the computation as name.SpectrCohTime.mat
%


function [y, f, t, phi] = SpectrCohTime(eeg, fMode, varargin)

[FreqRange, WinLengthSec, SR, Whiten, name] = ...
    DefaultArgs(varargin,{[0.1 50], 4, 1250, 1,[]});

switch fMode
    
    case 'c'
        
        
        WinLengthSample = 2^round(log2(WinLengthSec*SR));
        nFFT = 2*WinLengthSample;
        
        if Whiten ==1
            eeg = WhitenSignal(eeg,SR*2000,2);
        end
        
        tic
        [y, f, t, phi] = mtchglong(eeg,nFFT,SR,WinLengthSample,[],3,'linear',[],FreqRange);
        toc
        
        save([name '.' mfilename '.mat'], 'y','f','t','phi');
        
        
    case 'd'
        
        load([name '.' mfilename '.mat']);
        
        nChannels=size(y,3);
        
        figure('name', mfilename, 'NumberTitle','off');
        
        for Ch1=1:nChannels,
            for Ch2 = Ch1:nChannels
                subplot(nChannels, nChannels, Ch1 + (Ch2-1)*nChannels);
                if Ch1==Ch2
                    imagesc(t,f,20*log10(abs(y(:,:,Ch1,Ch2))+eps)');
                    axis xy; colormap(jet); %freezeColors; colorbar
                else
                    
                    imagesc(t,f,(abs(y(:,:,Ch1,Ch2)))');
                    axis xy; %coloration('dark'); freezeColors
                    % title(['Coherogram ' num2str(Ch1)]);
                    
                    
                    subplot(nChannels, nChannels, Ch2 + (Ch1-1)*nChannels);
                    
                    imagesc(t,f,squeeze(phi(:,:,Ch1,Ch2))');
                    axis xy; %coloration('phase'); freezeColors
                    % title(['Phasogram ' num2str(Ch1)]);
                end
            end
            % ForAllSubplots('set(gca, ''FontSize'', 6)')
        end
        xlabel('Time')
        ylabel('Frequency')
        
        
        
    case 'cd'
        
        
        WinLengthSample = 2^round(log2(WinLengthSec*SR));
        nFFT = 2*WinLengthSample;
        
        if Whiten ==1
            eeg = WhitenSignal(eeg,SR*2000,2);
        end
        
        tic
        [y, f, t, phi] = mtchglong(eeg, nFFT,SR,WinLengthSample,[],3,'linear',[], FreqRange);
        toc
        
        
        nChannels=size(y,3);
        
        figure('name', [mfilename ' - ' name], 'NumberTitle','off');
        
        for Ch1=1:nChannels,
            for Ch2 = Ch1:nChannels
                subplot(nChannels, nChannels, Ch1 + (Ch2-1)*nChannels);
                if Ch1==Ch2
                    imagesc(t,f,20*log10(abs(y(:,:,Ch1,Ch2))+eps)'); % Cross-spectrum
                    axis xy; colormap(jet);%freezeColors; colorbar
                else
                    
                    imagesc(t,f,(abs(y(:,:,Ch1,Ch2)))');  % Coherogram
                    axis xy; %coloration('dark'); freezeColors
                    
                    subplot(nChannels, nChannels, Ch2 + (Ch1-1)*nChannels);
                    
                    imagesc(t,f,squeeze(phi(:,:,Ch1,Ch2))'); % Phasogram
                    axis xy; %coloration('phase'); freezeColors
                    
                    %       #^$%&*(        if Ch2 == length(nChannels) & Ch1 == length(nChannels)
                    
                end
            end
            % ForAllSubplots('set(gca, ''FontSize'', 6)')
        end
        xlabel('Time')
        ylabel('Frequency')
        
    case 'cd_traces'
        
        
        WinLengthSample = 2^round(log2(WinLengthSec*SR));
        nFFT = 2*WinLengthSample;
        
        if Whiten ==1
            eeg = WhitenSignal(eeg,SR*2000,2);
        end
        
        tic
        [y, f, t, phi] = mtchglong(eeg, nFFT,SR,WinLengthSample,[],3,'linear',[], FreqRange);
        toc
        
        
        nChannels=size(y,3);
        
        figure('name', mfilename, 'NumberTitle','off');
        
        subplot 211
        plot(eeg); axis tight
        subplot 212
                           imagesc(t,f,20*log10(abs(y(:,:,1,1))+eps)'); % Cross-spectrum
                    axis xy; colormap(jet);%freezeColors; colorbar
              
                    
        xlabel('Time')
        ylabel('Frequency')

end
