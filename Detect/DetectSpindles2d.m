%function Spindle = DetectSpindles2d(FileBase, Channels,eSampleRate)
function Spindle = DetectSpindles2d(FileBase, varargin)

[Channels,eSampleRate] = DefaultArgs(varargin,{[1:32],1250});
% if isempty(GoodChannels) 
% 	GoodChannels=load([FileBase '.goodch']); 
% end
% GoodChannels=intersect(GoodChannels,NonSkippedChannels(FileBase));
Par = LoadPar([FileBase '.par']);
nT = FileLength([FileBase '.eeg'])/Par.nChannels/2;
nCh = length(Channels);

%parameters set:
neig =1;
win = round(eSampleRate*0.3); % window for SVD compression
Window = 100000;
MinISI = round(eSampleRate*0.5);
%NegPosLag = round(eSampleRate*0.3);
% % % ResampleTo = 250; %Hz
ResampleTo = eSampleRate; %Hz
FreqRange = [6 20]; % of spindles power
if mod(win,2)
	win = win+1;
	fprintf('making window size even\n');
end
nSeg = floor(nT/Window);
if mod(nT,Window)>0% adds one if there is one more chunk
	nSeg = nSeg +1;
end
Spindle = struct('t',[],'map',[]);
teigall = [];
for i=1:nSeg
	per = [(i-1)*Window+1:max(i*Window,nT)]';
	beg = (i-1)*Window+1;
	len = min(Window,length(per));
	seg = bload([FileBase '.eeg'], [Par.nChannels, len], beg*2*Par.nChannels);
	seg = seg(Channels,:)';% now seg is Time by channels
	seg = resample(seg,ResampleTo, eSampleRate);
	per = per(1:eSampleRate/ResampleTo:end);
    	%now filter
    	seg = ButFilter(seg,2,FreqRange/(ResampleTo/2),'bandpass');
	%now rectify and smooth the power
	nSmooth = round(ResampleTo/FreqRange(1)*2);
	if mod(nSmooth,2)==0
	nSmooth =nSmooth+1;
	end
	seg = Filter0(ones(nSmooth,1)/nSmooth, abs(seg));
	%FROM NOW ON IN ResampleTo sampling rate !!!!!!!!!!!!!!!!!!!!!!!
    	[teig, xeig] = PieceSVD(seg,neig,win);
	teig = sum(teig,2);
    	teigall = [teigall; teig];
    if 1
    %detect minima
    thr = 1*std(teig);
    pkspos = LocalMinima(-teig, MinISI, -thr);
    
    Spindle.t = [Spindle.t; per(pkspos(:))];
    Spindle.map = [Spindle.map; seg(pkspos,:)];
    end
 
    
end
thr = 1*std(teigall);
spdpos = LocalMinima(-teigall, MinISI, -thr)*eSampleRate/ResampleTo;

lab = [ones(length(Spindle.t),1); ones(length(spdpos),1)*2];
res = [Spindle.t; spdpos];

[res, ind ] = sort(res);
lab = lab(ind);
MakeEvtFile(res,[FileBase '.spd.evt'],lab,eSampleRate);

