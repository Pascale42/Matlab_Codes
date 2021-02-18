% function PowDP(FileBase, State, What, fMode,  :: FileName, FreqRange, resampling, SR, Channels, LogPlot)
%
% Computes Power Depth Profile of data
% FileBase :: FileBase of the dataset OR matrix containing the eeg data (=data mode) ONLY in 'compute mode'
% State :: THE or SWS
% What :: Name of the region/structure of which the PDP is computed
% fMode :: 'compute'  'display' or 'displaySpectra'
% DEFAULT VARIABLES ::
% FileName :: give the FileBase name if use the 'data mode' in FileBase otherwise leave it as []
% FreqRange :: [1 150]
% resampling :: to speed up computation is data is long enough, should be 2 or 10, use 1 for no resampling
% SR :: from Sampling Rate
% Channels :: channels to load if use 'FileBase mode'
% LogPlot :: plot PDP in log 1, otherwise 0



function PowDP(FileBase, State, What, fMode, varargin)

[FileName, FreqRange, resampling, SR, Channels, LogPlot] = ...
    DefaultArgs(varargin,{'',[1 150],1,1250,[],0});


switch fMode
    case 'compute'
        
        if ischar(FileBase)
            eeg=GetEegState(FileBase, Channels, 0, State);
            FileName=FileBase;
        else
            eeg=FileBase;
        end
        
        disp('...resample and whiten...')
        
        % Resample data if chosen
        reSampleRate=SR/resampling;
        eegr = resample(eeg,1,resampling);
        
        % Whitening
        weeg= WhitenSignal(eegr,reSampleRate*2000,1);
        
        % Parameters
        WinLengthSec=4;
        WinLengthSample = 2^round(log2(WinLengthSec*reSampleRate));
        nFFT = 2*WinLengthSample;
        
        disp('...computes power depth profile...')
        
        % Compute PDP for each channel
        [y,PDP.f]=mtptchd(weeg(:,1), [], [],nFFT,reSampleRate,WinLengthSample,[],3,'linear',[],FreqRange);
        PDP.y=zeros(length(PDP.f),size(weeg,2)); disp('ch 1 done');
        if size(weeg,2) > 1
            PDP.y(:,1)=y;
            for n=2:size(weeg,2); [y]=mtptchd(weeg(:,n),[],[],nFFT,reSampleRate,WinLengthSample,[],3,'linear',[],FreqRange);
                PDP.y(:,n)=y; disp(['ch ' num2str(n) ' done']);end
        else
            PDP.y(:,1)=y;
        end
        clear y n
        
        
        % Interpolating 50Hz and its harmonics
        FRange=floor(FreqRange(1) : FreqRange(2));
        FRange=FRange./50;
        a= FRange==floor(FRange);
        Fifties=FRange(a)*50;
        if Fifties(1)==0; Fifties=Fifties(2:end); end
        for n=1:length(Fifties)
            x=find(PDP.f >= Fifties(n)-1, 1, 'first');
            low=x-23; high=x+23;
            if low < 0; low=1; end
            if high > length(PDP.f); high=length(PDP.f); end
            y=interp1([PDP.f(low);PDP.f(high)],[PDP.y(low,:);PDP.y(high,:)],PDP.f(low:high),'pchip');
            PDP.y(low:high,:)=y;
            clear low high x
        end
        
        
        % For log scale
        if LogPlot
            PDP.Y=20*log10(PDP.y)+eps;
            PDP.zlog=imagesclog(PDP.f, 1:size(PDP.Y,2),PDP.Y,0);
        end
        
        
        
        if ischar(FileBase)
            save([FileBase '.PowDP.' State '.' What '.mat'],'PDP', 'FreqRange', 'resampling')
        else
            save([FileName '.PowDP.' State '.' What '.mat'],'PDP', 'FreqRange', 'resampling')
        end
        
    case 'display'
        
        if ~ischar(FileBase)
            FileBase=FileName;
        end
        load([FileBase '.PowDP.' State '.' What '.mat'])
        
        if LogPlot == 1
            figure('name', ['PDPs ' State '.' What],'NumberTitle','off');
            imagesc(PDP.f, 1:size(PDP.Y,2),PDP.zlog');colorbar;
            set(gca,'XScale','log'); xlabel('Frequency, Hz'); ylabel('channels'); colormap('jet')
            x=get(gca,'XTick'); set(gca,'XTickLabel',x); clear x; xlim(FreqRange)
        else
            figure('name', ['PDPs ' State '.' What],'NumberTitle','off');
            imagesc(PDP.f, 1:size(PDP.y,2),(20*log10(PDP.y(:,:))+eps)');
            set(gca,'yaxislocation','right'); set(gca,'YTick',[1:size(PDP.y,2)]) 
            colorbar; xlabel('Frequency, Hz'); ylabel('channels'); colormap('jet')
        end
        
    case 'displaySpectra'
        
        if ~ischar(FileBase)
            FileBase=FileName;
        end
        load([FileBase '.PowDP.' State '.' What '.mat'])
        
        for n=1:length(Channels);
            figure
            plot(PDP.f,PDP.Y(:,n)); grid on; axis tight
            ylabel('Power, dB'); xlabel('Frequency, Hz'); title(['Channel ' num2str(Channels(n))])
        end
        
        
    case 'displayPoly2'
        
        if ~ischar(FileBase)
            FileBase=FileName;
        end
        load([FileBase '.PowDP.' State '.' What '.mat'])
        Par=LoadPar([FileBase '.xml']);
        
        if LogPlot == 1
            figure('name', [FileBase ' - PDPs ' State '.' What],'NumberTitle','off');
            ch=1:size(PDP.Y,2);
            even = ch(2:2:length(ch));
            odd = ch(1:2:length(ch));
            subplot(1,2,1)
            imagesc(PDP.f, 1:length(even),PDP.zlog(:,even)');colorbar;
            set(gca,'XScale','log'); xlabel('Frequency, Hz'); ylabel('channels'); colormap('jet')
            x=get(gca,'XTick'); set(gca,'XTickLabel',x); clear x; xlim(FreqRange)
            subplot(1,2,2)
            imagesc(PDP.f, 1:length(odd),PDP.zlog(:,odd)');colorbar;
            set(gca,'XScale','log'); xlabel('Frequency, Hz'); ylabel('channels'); colormap('jet')
            x=get(gca,'XTick'); set(gca,'XTickLabel',x); clear x; xlim(FreqRange)
        else
            figure('name', [FileBase ' - PDPs ' State '.' What],'NumberTitle','off');
            ch=1:size(PDP.y,2);
            even = ch(2:2:length(ch));
            odd = ch(1:2:length(ch));
            subplot(1,2,1)
            imagesc(PDP.f, 1:length(even),(20*log10(PDP.y(:,even))+eps)');
            colorbar; xlabel('Frequency, Hz'); ylabel('channels'); colormap('jet')
            set(gca,'YTick',[1:length(even)]); set(gca,'YTickLabel',Par.AnatGrps(1).Channels(even)+1)
            subplot(1,2,2)
            imagesc(PDP.f, 1:length(odd),(20*log10(PDP.y(:,odd))+eps)');
            colorbar; xlabel('Frequency, Hz'); ylabel('channels'); colormap('jet')
            set(gca,'YTick',[1:length(odd)]); set(gca,'YTickLabel',Par.AnatGrps(1).Channels(odd)+1)
        end
%         savefig([FileBase '.PowDP.' State '.' What '.fig']);
%         close(gcf);

case 'display4shanks'
    
    if ~ischar(FileBase)
            FileBase=FileName;
        end
        load([FileBase '.PowDP.' State '.' What '.mat'])
        Par=LoadPar([FileBase '.xml']);
        
        figure('name', [FileBase ' - PDPs ' State '.' What],'NumberTitle','off');
            
            for n=1:4
            subplotfit(n, 4)
            imagesc(PDP.f, 1:length(Par.AnatGrps(n).Channels),(20*log10(PDP.y(:,8*n-7:8*n))+eps)');
            colorbar; xlabel('Frequency, Hz'); ylabel('channels'); colormap('jet')
            set(gca,'YTick',[1:length(Par.AnatGrps(n).Channels)]); set(gca,'YTickLabel',Par.AnatGrps(n).Channels+1)
            end
    

end















