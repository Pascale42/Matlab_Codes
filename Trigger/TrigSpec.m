%function [spec, f] = TrigSpec(X,T,Window, FreqRange,Fs,plotstyle)
% computes the triggered average power spectrum of 
% (so far) 1 channel signal X by time points T (sampling rate 
% should be the same as X) in +/- Window (in samples) from the trigger
% ouput spec =  log(TrigAvSpec)-log(WholeTraceAvSpec)
% so units are proportional to Db and reflect the difference
% it's easier to see the effect this way due to1/f shape of the spectrum
function [spec, f] = TrigSpec(X,T,varargin)
[Window, FreqRange, Fs,plotstyle] = DefaultArgs(varargin,{ ...
        2^10, [1 100], 1250, 'short'});

if length(Window)==1
    nSamples = 2*Window+1;
    Window = [Window Window];
else
    nSamples = sum(Window) + 1;
end
%Window = round(Window*Fs/1000);

%X=X(:); %deal with one dimension so far
X= WhitenSignal(X);
maxTime = size(X, 1);
nColumns = size(X,2);
MyTriggers = find(T>Window(1) & T<maxTime-Window(2));
BlockSize = floor(2000000/nSamples); % memory saving parameter
nTriggers = length(MyTriggers);
XTriggered = zeros(nTriggers,nSamples);
% go through triggers in groups of BlockSize to save memory
for Block = 1:ceil(nTriggers/BlockSize)
	BlockTriggers = MyTriggers(1+(Block-1)*BlockSize:min(Block*BlockSize,nTriggers));
	nBlockTriggers = length(BlockTriggers);
	
    TimeOffsets = repmat(-Window(1):Window(2), nBlockTriggers, 1);
	TimeCenters = repmat(T(BlockTriggers), 1, nSamples);
	TimeIndex = TimeOffsets + TimeCenters;
    Waves = X(TimeIndex);
	%Waves = reshape(Waves, [nBlockTriggers, nSamples, nColumns]);
	XTriggered(1+(Block-1)*BlockSize:min(Block*BlockSize,nTriggers),:) = Waves;
end
WrapXTrig = reshape(XTriggered',[1,nTriggers*nSamples]);

[yo, f]=mtcsglong(WrapXTrig,2^11,Fs,nSamples,0,2,'linear',[],FreqRange);
%yo=reshape(yo,[nTriggers,nSamples]);
TrigSpec = mean(yo,2);
[AllSpec, f] = mtcsd(X,nSamples,Fs,nSamples,0,5,'constant');
%keyboard
spec = log(TrigSpec) - log(AllSpec);
if (nargout ==0)
   % figure
   if (strcmp(plotstyle,'short'))
     plot(f,spec);set(gca,'XLim',[0 200]);
    else
        subplot(211)
        plot(f,log(TrigSpec),'g');
        hold on
        plot(f,log(AllSpec),'b');
        subplot(212)
        plot(f,spec);
    end
end
