% out={spec, f} = TrigStatistics(X,T,StatFun,Params,Window, Fs)
function out = TrigStatistics(X,T,StatFun,Params,varargin)

[Window, Fs] = DefaultArgs(varargin,{1000, 1250}); 
Window = nextpow2(Window*Fs/1000);

if length(Window)==1
    nWinSamples = 2*Window+1;
    Window = [Window Window];
else
    nWinSamples = sum(Window) + 1;
end
if size(X,1)<size(X,2)
    X=X'; 
end

nChannels = size(X,2);
nSamples = size(X,1);
MyTriggers = find(T>Window(1) & T<nSamples-Window(2));
nTriggers = length(MyTriggers);
XTriggered = zeros(nTriggers,nSamples);
TimeOffsets = repmat(-Window(1):Window(2), nTriggers, 1);
TimeCenters = repmat(T(MyTriggers), 1, nWinSamples);
TimeIndex = TimeOffsets + TimeCenters;
Waves = X(TimeIndex,:);
Waves = reshape(Waves', [nChannels,nTriggers,nWinSamples]);
Waves = permute(Waves, [2 3 1]);
out = feval('StatFun',Waves,Params);

