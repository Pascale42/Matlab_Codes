% function CohDP(data, State, :: fMode, FreqRange, resampling, Par, Channels, LogPlot, FileBase)
% 
%


function CohDP(data1, data2, State, varargin)

[fMode, FreqRange, resampling, Par, Channels, LogPlot, FileBase] = ...
    DefaultArgs(varargin,{'compute', [0.1 150],2,[],[],1,''});


switch fMode
    case 'compute'
        
%         if ischar(data)
%             eeg=GetEegState(data, Channels, 0, State);
%             Par = LoadPar([data '.xml']);
%         else
%             eeg=data; 
%         end
        
        % Resample data if chosen
        reSampleRate=Par.lfpSampleRate/resampling;
        eegr1 = resample(data1,1,resampling);
        eegr2 = resample(data2,1,resampling);
        
        % Whitening
        weeg1= WhitenSignal(eegr1,reSampleRate*2000,1);
        weeg2= WhitenSignal(eegr2,reSampleRate*2000,1);
        
        % Parameters
        WinLengthSec=4;
        WinLengthSample = 2^round(log2(WinLengthSec*reSampleRate));
        nFFT = 2*WinLengthSample;
        
        % Compute Coherence Depth Profile for each channel
        
        [y,CDP.f]=mtptchd([weeg1 weeg2(:,1)], [], [],nFFT,reSampleRate,WinLengthSample,[],3,'linear',[],FreqRange);
        CDP.y=zeros(length(CDP.f),size(weeg2,2));
        if size(weeg2,2) > 1
            CDP.y(:,1)=y;
            for n=2:size(weeg2,2); [y]=mtptchd([weeg1 weeg2(:,1)],[],[],nFFT,reSampleRate,WinLengthSample,[],3,'linear',[],FreqRange);CDP.y(:,n)=y(1,2); end
        else
            CDP.y=y(:,1,2); 
        end
        clear y n
        
        
        % Interpolating 50Hz and its harmonics
        FRange=floor(FreqRange(1) : FreqRange(2));
        FRange=FRange./50;
        a= FRange==floor(FRange);
        Fifties=FRange(a)*50;
        if Fifties(1)==0; Fifties=Fifties(2:end); end
        for n=1:length(Fifties)
            x=find(CDP.f >= Fifties(n), 1, 'first');
            low=x-23; high=x+23;
            if low < 0; low=1; end
            if high > length(CDP.f); high=length(CDP.f); end
            y=interp1([CDP.f(low);CDP.f(high)],[CDP.y(low,:);CDP.y(high,:)],CDP.f(low:high),'pchip');
            CDP.y(low:high,:)=y;
            clear low high x
        end
        
        
        % For log scale
        CDP.zlog=imagesclog(CDP.f, 1:size(CDP.Y,2),CDP.Y,0);
        
        
        
%          if ischar(data)
%             save([data '.PDP.' State '.mat'],'PDP', 'FreqRange', 'resampling')
%         else
            save([FileBase '.CohDP.' State '.mat'],'CDP', 'FreqRange', 'resampling')
%         end 
        
    case 'display'
        
        load([data '.CohDP.' State '.mat'])
        
        if LogPlot == 1
            figure('name', ['CohDP ' State],'NumberTitle','off');
            imagesc(CDP.f, 1:size(CDP.Y,2),CDP.zlog');colorbar;
            set(gca,'XScale','log'); xlabel('Frequency, Hz'); ylabel('channels')
            x=get(gca,'XTick'); set(gca,'XTickLabel',x); clear x; xlim(FreqRange)
        else
            figure('name', ['CohDP ' State],'NumberTitle','off');
            imagesc(CDP.f, 1:size(CDP.y,2),(20*log10(CDP.y(:,:))+eps)'); colorbar; xlabel('Frequency, Hz'); ylabel('channels')
        end
        
    case 'displaySpectra'
        
        load([data '.CohDP.' State '.mat'])
        
        for n=1:length(Channels);
            figure
            plot(CDP.f,CDP.Y(:,n)); grid on; axis tight
            ylabel('Coherence'); xlabel('Frequency, Hz'); title(['Channels ' num2str(Channels(n))])
        end
end















