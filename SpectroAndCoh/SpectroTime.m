% function SpectroTime(eeg, fMode, VARARGIN: FreqRange, WinLengthSec, SR, Whiten, name)
%
% fMode:  'c'  :: compute and save
%              'd'  :: display; only if saved before
%             'cd' ::compute and display, not saved
%
% DefaultArgs: FreqRange=[0.1 48]; WinLengthSec=4;
%              SR(sampling rate)=1250; Whiten=1; name(string)=[];
%
% Use 'name' to save the computation as name.SpectroTime.mat
%

function [y, f] = SpectroTime(eeg, fMode, varargin)

[FreqRange, WinLengthSec, SR, Whiten, name] = ...
    DefaultArgs(varargin,{[0.1 48], 4, 1250, 1, []});

switch fMode
    
    case 'c'
        
        WinLengthSample = 2^round(log2(WinLengthSec*SR));
        nFFT = 2*WinLengthSample;
        
        if Whiten ==1
            eeg = WhitenSignal(eeg,SR*2000,2);
        end
        
        %         tic
        %        [y, f, t] = mtchglong(eeg(:,1),nFFT,SR,WinLengthSample,[],3,'linear',[],FreqRange);
        %         toc
        
        %         if size(eeg,2) >1
        %             for n=2:size(eeg,2)
        tic
        %                  [Y, f, t] = mtchglong(eeg(:,n),nFFT,SR,WinLengthSample,[],3,'linear',[],FreqRange);
        [y, f, t] = mtchglong(eeg,nFFT,SR,WinLengthSample,[],3,'linear',[],FreqRange);
        toc
        %                 y=[y Y];clear Y
        %             end
        %         end
        
                   
        save([name '.' mfilename '.mat'], 'y','f','t');
        
        
        
    case 'd'
        
        load([name '.' mfilename '.mat']);
        nChannels=size(y,3);
        
        figure('name', [name ' - ' mfilename], 'NumberTitle','off');
        for n=1:nChannels
            subplot(nChannels,1,n)
            imagesc(t,f,20*log10(abs(y(:,:,n,n))+eps)');
                    axis xy; colormap(jet); colorbar
            ylabel('Frequency (Hz)'); xlabel('Time (s)');
            
            ForAllSubplots('set(gca, ''FontSize'', 8)')
        end
        
        
    case 'cd'
        
        WinLengthSample = 2^round(log2(WinLengthSec*SR));
        nFFT = 2*WinLengthSample;
        
        if Whiten ==1
            eeg = WhitenSignal(eeg,SR*2000,2);
        end
        
            
        tic
       [y, f, t] = mtchglong(eeg(:,1),nFFT,SR,WinLengthSample,[],3,'linear',[],FreqRange);
        toc
        
        if size(eeg,2) >1
            for n=2:size(eeg,2)
                tic
                 [Y, f, t] = mtchglong(eeg(:,n),nFFT,SR,WinLengthSample,[],3,'linear',[],FreqRange);
                toc
                y=[y Y];clear Y 
            end
        end
        
          nChannels=size(y,2);
        
        figure('name', [name ' - ' mfilename], 'NumberTitle','off');
        for n=1:nChannels
            subplot(nChannels,1,n)
            imagesc(t,f,20*log10(abs(y(:,n))+eps)');
                    axis xy; colormap(jet); %freezeColors; colorbar
            ylabel('Frequency (Hz)'); xlabel('Time (s)');
            
            ForAllSubplots('set(gca, ''FontSize'', 8)')
        end
        
        
end
