%function [DiffSpec, f, t, TrigSpec, BslnSpec] = ...
%           TrigSpecgram(X,T,Window,nWindows,Fs,Period, FreqRange, Whiten, StepFraction, NW)
% X - eeg trace, T  - times of the triggers
% Window - in samples, nWindos - number of bins, Fs - sample rate, 
% Period - periods where to sample the baseline spectra from
% FreqRange - [Fbeg Fend] in Hz
% Whiten =1 if you want to whiten the signal (removes 1/f slope)
% StepFraction - how much of the window you want to step in spectrogram
% calculation
function [DiffSpec, f, t, TrigSpec, BslnSpec] = TrigSpecgram(X,T,varargin)
[Window,nWindows,   Fs, Period, FreqRange, Whiten,StepFraction, NW] = DefaultArgs(varargin, ...
 { 2^10,       10,             1250, [],           [1 100] ,       1  ,        0.5 ,1 });

if length(Window)==1
    nSamples = (2*nWindows)*Window+1;
    WindowSize = Window;
    Window = [Window*nWindows Window*nWindows];
else
    %if the window is the vector
    nSamples = sum(Window);
    Window  = Window*nWindows;
end
if Whiten
    X = WhitenSignal(X);
end
if size(X,1)<size(X,2) X=X'; end
%X=X(:); %deal with one dimension so far
maxTime = size(X, 1);
nColumns = size(X,2);
MyTriggers = T(find(T>Window(1) & T<maxTime-Window(end)));
BlockSize = floor(2000000/nSamples); % memory saving parameter
nTriggers = length(MyTriggers);

if nTriggers==0
    error('no triggers found');
    spec =[];f=[];t=[];TrigSpec=[];AllSpec1=[];
    return
end


%WrapXTrig = reshape(XTriggered',[1,nTriggers*nSamples]);
nFFT = 2^10; %nSamples/nWindows;
swsX = SelectPeriods(X,Period,'c',0);
[BslnSpec, f] = mtpsd(swsX,nFFT,Fs,WindowSize,0,NW,'linear',[],FreqRange);
clear swsX;

TrigSpec =[];%zeros(length(f),nWindows*2-2,nColumns);
for i=1:nTriggers
    [yo, f,t]=mtpsg(X(MyTriggers(i)+[(-Window(1)+1):Window(2)],:),...
        nFFT,Fs,WindowSize,WindowSize-round(WindowSize*StepFraction), NW,'linear',[],FreqRange); %,NW,Detrend,nTapers);
    if (i==1)
        TrigSpec =yo;
    else
        TrigSpec = TrigSpec + yo;
    end
end

%t = t - t((nSamples-1)/2+1);

TrigSpec = TrigSpec/nTriggers;
TrigSpec = permute(TrigSpec,[2 1 3]);
% replicate average power spectra for all time bins
AllSpec = permute(repmat(BslnSpec,[1, 1, length(t)]),[3 1 2]);

%AllSpec = repmat(AllSpec1,[1, length(t) ,nColumns]);
%keyboard
DiffSpec = log(TrigSpec) - log(AllSpec);
DiffSpec = squeeze(DiffSpec);
AllSpec = squeeze(AllSpec);
TrigSpec = squeeze(TrigSpec);

t = t - median(t);
if (nargout ==0)
   % figure
   for i=1:nColumns
       subplot(nColumns,2,(i-1)*2+1);
       imagesc(t,f,squeeze(DiffSpec(:,:,i))');axis xy

       subplot(nColumns,2,(i-1)*2+2);
       imagesc(t,f,squeeze(TrigSpec(:,:,i))');axis xy

    end
end
