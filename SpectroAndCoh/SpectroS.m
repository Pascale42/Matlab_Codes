% function SpectroS(eeg, fMode, VARARGIN: FreqRange, WinLengthSec, SR, Whiten, name, tag)
%
% fMode:  'compute'  :: compute and save
%              'display'  ::  only if saved before
%              'cd' :: computes and displays; does not save
%
% DefaultArgs: FreqRange=[0.1 48]; WinLengthSec=4;
%              SR(sampling rate)=1250; Whiten=1; name(string)=[]; tag(name
%              channels)=[];
%
% Use 'name' to save the computation as name.SpectroTime.mat
%

function [y, f] = SpectroS(eeg, fMode, varargin)

[FreqRange, WinLengthSec, SR, Whiten, name, tag] = ...
    DefaultArgs(varargin,{[0.1 48], 4, 1250, 1, [], []});

switch fMode
    
    case 'compute'
        
        WinLengthSample = 2^round(log2(WinLengthSec*SR));
        nFFT = 2*WinLengthSample;
        
        if Whiten ==1
            eeg = WhitenSignal(eeg,SR*2000,2);
        end
        
        % time
        tic; [y, f, t] = mtchglong(eeg,nFFT,SR,WinLengthSample,[],3,'linear',[],FreqRange); toc
        
        % Freq
        tic;   [Y, F]= mtchd(eeg(:,1),nFFT,SR,WinLengthSample,[],3,'linear',[],FreqRange);
        if size(eeg,2) >1
            for n=2:size(eeg,2)
                tic
                [Y1]= mtptchd(eeg(:,n),[],[],nFFT,SR,WinLengthSample,[],3,'linear',[],FreqRange);
                toc
                Y=[Y Y1];clear Y1
            end
        end
        toc
        
        
        save([name '.' mfilename '.mat'], 'y','f','t','Y','F');
        
        
        
    case 'display'
        
        load([name '.' mfilename '.mat']);
        nChannels=size(y,3);
        
        figure('name', [name ' - ' mfilename ' - ' tag], 'NumberTitle','off');
        for n=1:nChannels
            % Time
            subplot(nChannels,4,[(4*n-3): (4*n-1)])
            imagesc(t,f,20*log10(abs(y(:,:,n,n))+eps)');
            axis xy; colormap(jet); colorbar
            ylabel('Frequency (Hz)'); xlabel('Time (s)');
            
            % Freq
            subplot(nChannels,4,4*n)
            plot(F,20*log10(abs(Y(:,n))+eps),'r'); grid on; % power spectrum density
            ylabel('psd (dB)'); xlabel('Frequency (Hz)');
            
            ForAllSubplots('set(gca, ''FontSize'', 8)')
        end
        
    case 'cd'
        
        WinLengthSample = 2^round(log2(WinLengthSec*SR));
        nFFT = 2*WinLengthSample;
        
        if Whiten ==1
            eeg = WhitenSignal(eeg,SR*2000,2);
        end
        
        % time
        tic; [y, f, t] = mtchglong(eeg,nFFT,SR,WinLengthSample,[],3,'linear',[],FreqRange); toc
        
        % Freq
        tic;   [Y, F]= mtchd(eeg(:,1),nFFT,SR,WinLengthSample,[],3,'linear',[],FreqRange);
        if size(eeg,2) >1
            for n=2:size(eeg,2)
                tic
                [Y1]= mtptchd(eeg(:,n),[],[],nFFT,SR,WinLengthSample,[],3,'linear',[],FreqRange);
                toc
                Y=[Y Y1];clear Y1
            end
        end
        toc
        
        nChannels=size(y,3);
        figure('name', [name ' - ' mfilename ' - ' tag], 'NumberTitle','off');
        for n=1:nChannels
            % Time
            subplot(nChannels,4,[(4*n-3): (4*n-1)])
            imagesc(t,f,20*log10(abs(y(:,:,n,n))+eps)');
            axis xy; colormap(jet); colorbar
            ylabel('Frequency (Hz)'); xlabel('Time (s)');
            
            % Freq
            subplot(nChannels,4,4*n)
            plot(F,20*log10(abs(Y(:,n))+eps),'r'); grid on; % power spectrum density
            ylabel('psd (dB)'); xlabel('Frequency (Hz)');
            
            ForAllSubplots('set(gca, ''FontSize'', 8)')
        end
        
end
