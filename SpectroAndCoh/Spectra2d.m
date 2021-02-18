%function Spec = Spectra2d(FileBase, Window,FreqRange,GoodChannels,eSampleRate)
% WIndow in secs
function Spec = Spectra2d(FileBase, varargin)

[Window,FreqRange,GoodChannels,eSampleRate] = DefaultArgs(varargin,{0.5,[1 100],[1:32],1250});
if isempty(GoodChannels) 
	GoodChannels=load([FileBase '.goodch']); 
end
GoodChannels=intersect(GoodChannels,NonSkippedChannels(FileBase));
Par = LoadPar([FileBase '.par']);
nT = FileLength([FileBase '.eeg'])/Par.nChannels/2;
nCh = length(GoodChannels);
ResampleTo = 250;
BufWindow = 100000;	
WindowSamples = 2^round(log2(Window*ResampleTo)+1);
FreqRes = 0.5;
nFFT = 2^round(log2(ResampleTo/2/FreqRes));
Spec = struct('y',[],'t',[],'f',[]);

nSeg = floor(nT/BufWindow);
if mod(nT,BufWindow)>0% adds one if there is one more chunk
	nSeg = nSeg +1;
end
h = waitbar(0,'Please wait...');
maxt = 0;
for i=1:nSeg
	waitbar(i/nSeg,h);
	per = [(i-1)*BufWindow+1:max(i*BufWindow,nT)]';
	beg = (i-1)*BufWindow+1;
	len = min(BufWindow,length(per));
	seg = bload([FileBase '.eeg'], [Par.nChannels, len], beg*2*Par.nChannels);
	seg = seg(GoodChannels,:)';% now seg is Time by channels
	seg = resample(seg,ResampleTo, eSampleRate);
	seg = WhitenSignal(seg,[],1);
	[y,f,t]= mtcsglong(seg,nFFT,ResampleTo,WindowSamples,[],2.5,'linear',[],FreqRange);
	dt = t(2)-t(1);
	Spec.y = cat(1,Spec.y,y);
	Spec.t = [Spec.t(:) ;t(:)+dt/2+maxt];
	maxt = Spec.t(end);
end
Spec.f = f;
close(h);


	
	
