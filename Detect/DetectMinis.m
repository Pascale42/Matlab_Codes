%function Peaks =  DetectMinis(eeg, FreqRange, Sign, Thresh, SampleRate, MinInter, Periods, Verbose)
% detects minis from eeg of Sign (+/- 1) in FreqRange at Tresh that are more than MinInter msec apart 
% verbose goes through chunks and plots the detected ones.
% Periods restrict detection to particular intervals
function Peaks = DetectMinis(eeg, varargin)
[FreqRange, Sign, Thresh,  SampleRate, MinInter, Periods , Verbose] = DefaultArgs(varargin,  ...
    {[50 200],        1,       2,             1250,            10,       [] ,                0     });
MinInter = ceil(MinInter*SampleRate/1000); % minial interval for the detection in samples
BuffSize = 10; %seconds of display buffer
nClu = 5; % number of clusters to split mini's
RunMeanWin=1; % sec - window for running mean computation
% b = fir1(20, 10/SampleRate*2,'high');
% seeg =Filter0(b, eeg);
win = RunMeanWin*SampleRate;
if ~mod(win,2) win =win+1; end
seeg = eeg - RunningMean(eeg,win);

isok=0;                                                                     %DETECTION BEGINS HERE
while ~isok
    
    % first filter the trace
    [b,a] = butter(4,FreqRange/SampleRate*2);
    feeg = filtfilt(b,a,eeg);

    FiltPeaks = SignalPeaks(feeg, Sign, Thresh, MinInter, 'log', Periods,0);
    if 0
        figure
        PeakAmp =seeg(FiltPeaks);
        if min(PeakAmp)<0
            shift = abs(min(PeakAmp))+1;
            PeakAmp = PeakAmp + shift ;
        else 
            shift =0;
        end
        hist(log10(Sign*PeakAmp),100);
        keyboard
    end
    EegPeaks =  SignalPeaks(seeg, Sign, Thresh*0.8, MinInter, 'log', Periods,0);
    if isempty(EegPeaks) | isempty(FiltPeaks)
        Peaks = [];
        return;
    end
    [EegPeaksInd, FiltPeaksInd ] = CloserThen(EegPeaks, FiltPeaks, fix(MinInter));
    Peaks = FiltPeaks(FiltPeaksInd);
    fprintf('Detected %d events\n',length(Peaks));
    Amp = feeg(Peaks);
    lAmp = log(abs(Amp));
    %sort mini by amplitude
    [dummy AmpClu] = histcI(lAmp,linspace(min(lAmp), max(lAmp),nClu));
    
    if Verbose
        figure(777)
        Buff = SampleRate*BuffSize;
        nTime = length(eeg);
        nSeg = fix(nTime/Buff)+1;
        for i=1:nSeg
            Seg = [1+Buff*(i-1) : min(Buff*i, nTime)];
            Seg = Seg(:);
            EegSeg = seeg(Seg);
            myPeaks = Peaks(find(Peaks>Seg(1) & Peaks<Seg(end)));
            if isempty(myPeaks)
                continue;
            end
            plot(Seg*1000/SampleRate, EegSeg);
            hold on
            Lines(myPeaks*1000/SampleRate,[],AmpClu);
            hold off
            pause
        end
    end
    if Verbose
        isok = ~input('do you want to redetect (0 or 1)?\n');
        if ~isok & Verbose
            fprintf('current Thresh = %2.2f, MinInter = %d, FreqRange = [%d %d]\n', Thresh, MinInter,FreqRange);
            tinp = input('enter new: Thresh   MinInter  FreqRange :  ','s');
            tinp = str2num(tinp);
            Thresh = tinp(1); MinInter = tinp(2); FreqRange = tinp(3:4);
        end
    else
        isok = 1;
    end
end

Peaks  = [Peaks(:) Amp(:) AmpClu(:)];
