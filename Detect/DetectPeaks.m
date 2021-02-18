%function Out = DetectPeaks(FileBase,PeakSign,MinCycles,Thresh, visualize, 
%                       SampleRate,Periods,FilterType,FreqRange,Chans, FileExt, Control)
% detects the peaks of oscillations from the power and filtered traces computed
% by PowerComputeInFile 
% Input:  every entry could be a vector or one value for all
%           PeaksSign - +1 or -1 for peaks/troughs
%           MinCycles - min number of cycles to include in analysis 
%           Thresh - log power threshold in number of std
%           Periods - allows to do detection just from that region
%	    FileExt - (def =.eeg) of teh file to use
% Output: structure array for each signal (chan/freq.range): 
% Out(i). 
%           Bursts = [BurstCenterTime, BurstLength, BurstPow, BurstFreq BurstBeg BurstEnd]
%           Peaks = [PeakTime, PeakPosInBurst, PeakPow, FromWhichBurst, OutOf]
%           BurstQuality = [psd2isi psd2isi2 psdpeak psdsnr AvPeaksSnr]; - look at code for meaning
%           nBursts, nPeaks, FreqRange, Chan, Thresh, MinCycles 

function Out = DetectPeaks(FileBase,varargin)

% this is the minmal wave amplitute that will be taken for peak detection
MINWAVETHR = 0.1; %number of standard deviations for wave detection
MINPOW = 0; % minimal power to consider in analysis
FILTORDER = 400; % order for non wavelet filter

% parameters for CleanPeaks
MAXASSYMETRY = 4; %  maximal assimetry of the peak slopes
HALFAMP = 0.3; % at what part of the peak measure the derivatives
MINPEAKSNR = 0.5; % SNR for the peak amp (from 0), to 
PsdSnrThr = 0;

% set the default input parameters
[PeakSign,  MinCycles,Thresh, visualize, SampleRate, Periods,FilterType,FreqRange,Chans, FileExt, Control ] = DefaultArgs(varargin, ...
    {     -1,            2,           1,                 0,          1250,             [],     'butt',         [],               [],        'eeg',      1       });

if ~isempty(FreqRange) & ~isempty(Chans)
    Chans = GetChannels(FileBase, Chans,0);
    FreqRange = GetFreqRange(FreqRange,length(Chans));
end
% load Chans and Freq info or make new
oldFreqRange=[];oldChans=[];
if strcmp(FilterType,'wav') & FileExists([FileBase '.osc.mat'])
    if ~isempty(FreqRange) & ~isempty(Chans)
        oldFreqRange = FreqRange; oldChans = Chans;
    end
    load([FileBase '.osc.mat']);
    if ~isempty(oldFreqRange) & ~isempty(oldChans)
        [oldFreqRange, FreqRange] = swap(oldFreqRange, FreqRange);
        [oldChans, Chans] = swap(oldChans, Chans);
    end
elseif isempty(FreqRange) & isempty(Chans)
    fprintf('need to know FreqRange and Chans for that setting\n');
    return;
end

nCh = length(Chans);
nSignals = sum(cellfun('prodofsize',FreqRange)/2);
par = LoadPar([FileBase '.par']);

if length(Thresh)==1   Thresh = Thresh*ones(nSignals,1); end
if length(PeakSign)==1   PeakSign = PeakSign*ones(nSignals,1); end
if length(MinCycles)==1   MinCycles = MinCycles*ones(nSignals,1); end

% if strcmp(FileExt,'csd')
%     load([FileBase '.csd.ch'],'-MAT');
%     nChannels = length(CsdChannels{1});
% else
    nChannels = par.nChannels;
%end

for ch=1:nCh
    myCh = Chans(ch);
    nFreqRange = size(FreqRange{ch},1);
    eeg = readsinglech([FileBase '.' FileExt], nChannels,myCh);
    for freq =1:nFreqRange
        myFreqRange = FreqRange{ch}(freq,:);
        sig = (ch-1)*nFreqRange+freq; % index of the signal
        fprintf('Detecting oscillations in [%d %d] Hz range on channel %d\n',myFreqRange,myCh);
        
        if strcmp(FilterType, 'wav' ) & FileExists([FileBase '.oscpow'])
            pow = readsinglech([FileBase '.oscpow'], nSignals, sig);
            wav = readsinglech([FileBase '.oscwav'], nSignals, sig);
        else
            if strcmp(FilterType,'fir1')
                FirFiltb = fir1(FILTORDER,myFreqRange/SampleRate*2);
                wav = Filter0(FirFiltb, eeg);
		pow = abs(wav);
		FirFiltb = fir1(FILTORDER,myFreqRange.*[1 1.5]/SampleRate*2);
                wav = Filter0(FirFiltb, eeg);
            elseif strcmp(FilterType,'butt')
                [b,a] = butter(2,myFreqRange/SampleRate*2);
                wav = filtfilt(b,a,eeg);
		pow = abs(wav);
		[b,a] = butter(2,myFreqRange.*[1 1.5]/SampleRate*2);
                wav = filtfilt(b,a,eeg);
            else
                wav = filtereeg(eeg,myFreqRange,[],SampleRate,[],10);
		pow = abs(wav);
                %wav = eegfilt(eeg,SampleRate,FreqRange(1),FreqRange(2));
            end
            
	    
%              % smooth the ripple
%              FiltLength = round(SampleRate/myFreqRange(1));
%              FiltLength = FiltLength+mod(FiltLength,2); % make it odd
%              pow = Filter0(ones(FiltLength,1)/FiltLength,pow);
%              
        end
        % now smooth the power 
        % smooth with the square window of 1/4 the cycle length
        AvgFiltOrder = round(SampleRate/(myFreqRange(1)));
        AvgFiltOrder=AvgFiltOrder+(mod(AvgFiltOrder,2)-1); 
        AvgFiltb = ones(AvgFiltOrder,1)/AvgFiltOrder;
        pow = Filter0(AvgFiltb, pow);
        
        %find times of nozero power within Periods
        GoodPow = find(pow>MINPOW);
        if ~isempty(Periods)
            %        GoodPow = WithinRanges(GoodPow, Periods);
            GoodPow = SelectPeriods(GoodPow,Periods,'d',1);
        end
        
        % compute the threshold from log transformed power that is greater then MINPOW
        AvPow = mean(log(pow(GoodPow)));
        StdPow = std(log(pow(GoodPow)));
        
        myThresh = Thresh(sig);
        myMinCycles = MinCycles(sig);
        % now do the detection untill the user is satisfied %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        isok=0;                                                                     %REDETECTION BEGINS HERE
        while ~isok
            
            PowThr = myThresh*StdPow;
            PowThr =AvPow + PowThr;
            
            % may do some further cut of the low power
            %        p(find(p<-3))=0;
            %         p = (p-mean(p))/std(p);
            
            % computes the actual power peak freq, of the filtered trace
           % [pxx, pcx, f] =Spsd(wav,2^9,SampleRate);
           % [m mi ] =max(pxx);
          %  PeakFreq = f(mi);
           PeakFreq = mean(myFreqRange);
	    
            % detects peaks of a given sign in filtered signal , that
            % are apart by more than 1/2 periods of selected frequency (may not need this)
            % and less than MINWAVETHR std of the filtered signal 
            WavThr = -PeakSign(sig)*std(wav(GoodPow))*MINWAVETHR; %- mean(wav(GoodPow));
            MeanPeakInt = SampleRate/PeakFreq;
            BetweenPeaksInt = MeanPeakInt/2;
            
            %get peaks within GoodPow indexing
            Peaks = LocalMinima(-PeakSign(sig)*wav(GoodPow), BetweenPeaksInt, WavThr );
            
            %get peak times within GoodPow in the total indexing
            GoodPeaks = GoodPow(Peaks);
%commented            
            %get those troughs that have bigger power than threshold
%           HighPeaks = GoodPeaks( find( log(pow(GoodPeaks)) > PowThr ) );
%           HighPeaks = HighPeaks(:);
            HighPeaks = Peaks(:);
            % now analyze the duration of peaks bursts etc 
            % quick but not safe as the freq. may change for different bursts
            %        [Burst, BurstLen, SpkPos, OutOf] = SplitIntoBursts(HighPeaks, 3*MeanPeakInt);
            
            % laborous sorting of peaks bursts:
            % get the borders of the power bursts just below the threshold for the peaks
            %BurstBoundary = ThreshCross( pow, exp(PowThr*0.8), MeanPeakInt);  
	    
	    MinBurstDuration = SampleRate*myMinCycles/myFreqRange(2);
            BurstBoundary = ThreshCross( pow, exp(PowThr), MinBurstDuration);  
	    %keyboard
	    %            BurstBoundary = intersect(GoodPow,BurstBoundary(:));
            %           BurstBoundary = reshape(BurstBoundary,2,[])';
            
            %test block
            if 0
                eegplot(unity([eeg(:) pow(:) wav(:)])',1250);
                keyboard
            end
            nBurst = size(BurstBoundary,1);
            OutPeaks = []; OutBursts = []; OutBurstQuality=[];
            BurstCnt =1;
            
	    % now go through bursts and sort them out :))
            for i=1:nBurst
                % get peaks withing this boundary
                myPeaks = HighPeaks(find(HighPeaks>BurstBoundary(i,1) & HighPeaks<BurstBoundary(i,2)));
                seg = [BurstBoundary(i,1):BurstBoundary(i,2)];
                
                eegseg = eeg(seg);
                powseg = pow(seg);
                wavseg = wav(seg);
                [b,a] = butter(2,100/SampleRate*2);
                aveegseg = filtfilt(b,a,eegseg);
%                aveegseg = Filter0(AvgFiltOrder, eegseg);
                % check each peak , get rid of bad ones
                % cleaning paramters CleanPars are taken from defaults, but can be changed 
                CleanPars = [MAXASSYMETRY, HALFAMP, MINPEAKSNR];
                
                [myCleanPeaks PeaksSnr] = CleanPeaks(myPeaks-BurstBoundary(i,1) +1, wavseg, aveegseg, PeakSign(sig), CleanPars);
%                myCleanPeaks =  myPeaks-BurstBoundary(i,1) +1;
		%PeaksSnr = 0;
		myPeaks = BurstBoundary(i,1) -1 + myCleanPeaks;
                
                
                if length(myPeaks)< myMinCycles
                    continue;
                else
                    AvPeaksSnr = mean(PeaksSnr);
                    % find the central peak - one that has biggest power
                    [dumb centind] = max(pow(myPeaks));
                    myBurstCentPeak = myPeaks(centind);
                    myBurstPow = mean(log(powseg));
                    myBurstLen = length(myPeaks);
                    myBurstFreq = SampleRate/median(diff(myPeaks)); % in Hz
                    myPeaksPos = [1: myBurstLen] - centind;
                    myPeaksPow = log(pow(myPeaks));
                    [psdseg, f]=mtcsd(eegseg,2^10,SampleRate,2^fix(log2(length(seg)/2)),[],1.5,'linear',[], [0.8 1.2].*myFreqRange(:)');
                    psdpeaks = LocalMinima(-log(psdseg), -mean(psdseg)+2*std(psdseg), length(find(f<6)));
                    
                    %find psd peaks within our freq. of interest
                    %myPsdPeak = WithinRanges(f(psdpeaks),myFreqRange.*[0.8 1.2]); 
                    myPsdPeak = f(psdpeaks);
                    myPsdPeak = psdpeaks(find(myPsdPeak));
                    [psdmax pind] = max(psdseg(myPsdPeak));
                    myPsdPeak = f(myPsdPeak(pind)); % finally have the freq. of the peak 
                    psd2isi = abs(1 - myPsdPeak / myBurstFreq);
                    psdpeak =psdmax / mean(psdseg);
                    % compute the snr for the psd
                    fin = find(f>myFreqRange(1) & f<myFreqRange(2)); fout = setdiff([1:length(f)],fin);
                    psdsig1 = std(psdseg(fin));
                    psdsig0 = std(psdseg(fout));
                    psdsnr = 10*log10(psdsig1/psdsig0);
                    
                    if isempty(psd2isi)
                        psd2isi=-1;
                        psdpeak = 0;
                    end
                    
                    [dumb BigPsdPeak] = max(psdseg);
                    BigPsdPeak = f(BigPsdPeak); %  freq. of the biggest peak
                    psd2isi2 = abs(1 - BigPsdPeak / myBurstFreq);
                    Duration = diff(BurstBoundary(i,:),1,2);             
                    
                    %one more selection: by PsdSnr
                    if psdsnr < PsdSnrThr 
                        continue;
                    end
                    
                    %test block 
                    if visualize
                        fprintf('Burst #%d : power %2.2f; psd2isi %2.2f; psdsnr %2.2f; psdpeak %2.2f\n', ...
                            BurstCnt, (myBurstPow-AvPow)./ StdPow, psd2isi,psdsnr, psdpeak);
                        figure(666)
                        clf
                        subplot(3,1,[1 2]); % plots traces and peaks
                        Append = fix(0.1*Duration);
                        seg = [max(1,BurstBoundary(i,1)-Append):min(BurstBoundary(i,2)+Append, length(pow))];
                        powseg = pow(seg);
                        wavseg = wav(seg);
                        eegseg = eeg(seg);
                        plot(seg*1000/SampleRate,unity([wavseg(:) eegseg(:) powseg(:)]));
                        ax = axis;
                        hold on
                        Lines(myPeaks*1000/SampleRate, ax(3:4),'r');
                        hold off
                        
                        subplot(3,3,7);
                        plot([2:myBurstLen],diff(myPeaks)*1000./SampleRate,'+');
                        ylim(MeanPeakInt*[0.2 3]*1000/SampleRate);
                        xlim([1 myBurstLen+1]);
                        
                        subplot(3,3,8);
                        plot([1:myBurstLen],wav(myPeaks),'+');
                        xlim([0 myBurstLen+1]);
                        
                        subplot(3,3,9);
                        plot(f, 10*log(psdseg));
                        ax = axis;
                        hold on
                        Lines(myPsdPeak,ax(3:4),'r');
                        hold off
                        %xlim(PeakFreq*[0.3 3]);
                        if Control                        
                            pause
                        end
                    end
                    OutPeaks = [OutPeaks; myPeaks(:) myPeaksPos(:) myPeaksPow(:) ...
                            ones(myBurstLen,1)*BurstCnt ones(myBurstLen,1)*myBurstLen ];
                    
                    OutBursts(BurstCnt,:) = [ myBurstCentPeak myBurstLen myBurstPow myBurstFreq BurstBoundary(i,:) myPsdPeak];
                    
                    OutBurstQuality(BurstCnt, :) = [psd2isi psd2isi2 psdpeak psdsnr AvPeaksSnr];
                    
                    BurstCnt = BurstCnt+1;
                end % end of if visualize
            end % end of for BurstsBoundary loop
            % normalize the pow to the z-score
            if BurstCnt>1
                OutPeaks(:,3) = (OutPeaks(:,3) - AvPow)./ StdPow;
                OutBursts(:,3) = (OutBursts(:,3) - AvPow)./ StdPow;
            end
            if Control
                if size(OutBursts,1)==0
                    break;
                else
                    figure(666)
                    Labels = {'psd2isi','psd2isi2', 'psdpeak', 'psdsnr', 'AvPeaksSnr'; ...
                            'BurstLength(cycles)', 'BurstPower (std)', 'BurstFrequencyISI(Hz)','BurstFrequencySpec(Hz)', 'BurstDuration(msec)'}; 
                    HistMatrix(permute(cat(3,OutBurstQuality(:,[1:5]),...
                        [OutBursts(:,[2:4 7]), diff(OutBursts(:,5:6),1,2)*1000/SampleRate]),[1 3 2]), max(10,fix(BurstCnt/5)), [],Labels);
                    
                    isok = ~input('do you want to redetect (0 or 1)?\n');
                    if ~isok
                        fprintf('current Thresh = %2.2f, MinCycles = %d, PsdSnrThr = %2.2f\n', myThresh, myMinCycles, PsdSnrThr);
                        tinp = input('enter new: Thresh  MinCycles PsdSnrThresh :  ','s');
                        tinp = str2num(tinp);
                        myThresh = tinp(1); myMinCycles = tinp(2); PsdSnrThr = tinp(3);
                    end
                end
            else
                isok=1;
            end
            % and show goes on again if you have not accepted what I detected :))
        end       % end of while ~isok
        
        Out(sig).Bursts = OutBursts;
        Out(sig).Peaks = OutPeaks;
        Out(sig).BurstQuality = OutBurstQuality;
        Out(sig).nBursts = size(OutBursts,1); 
        Out(sig).nPeaks = size(OutPeaks,1);         
        Out(sig).FreqRange = myFreqRange;
        Out(sig).Chan = myCh;
        Out(sig).Thresh = myThresh;
        Out(sig).MinCycles = myMinCycles;
    end
end

if nargout<1
    save([FileBase '.DetectPeaks.mat'], 'Out');
end