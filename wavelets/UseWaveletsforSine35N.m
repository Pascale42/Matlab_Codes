function [proDomFreq,proMeanPower,StimFreqVect] = UseWaveletsforSine35N(FileBase,GoodChannels,WofI,varargin)

[Period,FreqRange,normalize] = DefaultArgs(varargin,{[],[6 13]},0);

% calculates dominant freq, mean power from wavelet
% transform
% WofI: window of interest, in secs [time before stim, duration of stim,
% time after stim)
% FreqRange, default 6-12Hz

% averages the results over the number of stim repetitions



Par = LoadPar([FileBase '.xml']);
Eeg = LoadBinary([FileBase '.eeg'],1,Par.nChannels);

if isempty(Period)
    Period = [1 length(Eeg)];
elseif ischar(Period)
    Period = [1 ceil(length(Eeg)./ Par.lfpSampleRate)];
end
clear Eeg


%PerMat = Period * Par.lfpSampleRate;
StimOnBegin = load([FileBase '.' num2str(Period(1)) '-' num2str(Period(2)) '.StB.evt']);
StimFreqVect = StimOnBegin(:,2);
StimOnBegin = StimOnBegin(:,1) *Par.lfpSampleRate / 1000;   %gets rid of tag and turns evts time (saved in ms in Evtfile) into datapoints

%disp(['n=' num2str(length(StimOnBegin)) 'stim repeats'])    % so we have the n corresponding to the graph

WofI = WofI .* Par.lfpSampleRate;

proDomFreq = zeros(length(StimOnBegin),(sum(WofI)+1),length(GoodChannels));
proMeanPower = zeros(length(StimOnBegin),(sum(WofI)+1),length(GoodChannels));


for i=1:length(GoodChannels)
    disp(['Channel # ' num2str(GoodChannels(i)) ' , out of ' num2str(length(GoodChannels))])
    waveletfilename = ['Wavelet' FileBase '_Ch' num2str(GoodChannels(i)) '.mat'];
    if file_exists(waveletfilename)
        disp('loading wavelets...')
        load(waveletfilename);
    else
        disp('calculating wavelets...')
        [wave,period,scale,coi] = makeBigWaveletMatrix(FileBase,GoodChannels(i));  %appelle la fonction pour faire les wavelettes
    end
    
    WinFreq = find(scale < (1./FreqRange(1)) & scale > (1./FreqRange(2)));
    
    
    for j=1:length(StimOnBegin)
        disp(['stim # ' num2str(j) ' , out of ' num2str(length(StimOnBegin))])
        tempWofI = [StimOnBegin(j)-WofI(1) StimOnBegin(j)+WofI(2)+WofI(3)];
        power = abs(wave(:,tempWofI(1):tempWofI(2))) .* abs(wave(:,tempWofI(1):tempWofI(2)));
        [~,preproDomFreq] = max(power,[],1);
        proDomFreq(j,:,i) = 1./scale(preproDomFreq)';
        
        proMeanPower(j,:,i) = sum(power(WinFreq(1):WinFreq(end),:),1);
        
    end
        
    
end
        
        
        
        
        
