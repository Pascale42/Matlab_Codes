%function Delta = DetectDeltaSVD(FileBase, GoodChannels,eSampleRate)
function Delta = DetectDeltaSVD(FileBase, varargin)
Par = LoadPar([FileBase '.xml']);
[GoodChannels,eSampleRate] = DefaultArgs(varargin,{[1:Par.nChannels],1250});
if isempty(GoodChannels) 
	%GoodChannels=load([FileBase '.goodch']); 
 GoodChannels =  NonSkippedChannels(FileBase);
end
%GoodChannels=intersect(GoodChannels,NonSkippedChannels(FileBase));

nT = FileLength([FileBase '.eeg'])/Par.nChannels/2;
nCh = length(GoodChannels);

%parameters set:
neig =1;
win = round(eSampleRate*0.3); % window for SVD compression
Window = 100000;
MinISI = round(eSampleRate*0.1);
NegPosLag = round(eSampleRate*0.3);
ResampleTo = 25; %Hz

if mod(win,2)
	win = win+1;
	fprintf('making window size even\n');
end
nSeg = floor(nT/Window);
if mod(nT,Window)>0% adds one if there is one more chunk
	nSeg = nSeg +1;
end
Delta = struct('neg',struct('t',[],'map',[]),'pos',struct('t',[],'map',[]));
teigall = [];
for i=1:nSeg
	per = [(i-1)*Window+1:max(i*Window,nT)]';
	beg = (i-1)*Window+1;
	len = min(Window,length(per));
	seg = bload([FileBase '.eeg'], [Par.nChannels, len], beg*2*Par.nChannels);
	seg = seg(GoodChannels,:)';% now seg is Time by channels
	seg = resample(seg,ResampleTo, eSampleRate);
	per = per(1:eSampleRate/ResampleTo:end);
    	%now filter
    	seg = ButFilter(seg,2,10/(ResampleTo/2),'low');
	%FROM NOW ON IN ResampleTo sampling rate !!!!!!!!!!!!!!!!!!!!!!!
    	[teig, xeig] = PieceSVD(seg,neig,win);
	teig = sum(teig,2);
    	teigall = [teigall; teig];
    if 1
    %detect minima
    thr = 1*std(teig);
    pkspos = LocalMinima(-teig, MinISI, -thr);
    pksneg = LocalMinima(teig, MinISI, -thr);
	    
%      pksall = [pksneg(:); pkspos(:)];
%      cluall = [ones(length(pksneg),1); ones(length(pkspos),1)*2];
%      [pksall, si ] = sort(pksall);
%      cluall= cluall(si);
%      
%      [ccg,tbin,pairs] = CCG(pksall,cluall,NegPosLag,0,ResampleTo,[1 2], 'scale');
%      
%      % now in pairs(:,1) is the index for which pksneg has pkspos within NegPosLag 
%      
%      bothi = intersect(pairs(:,1), find(cluall==1));
        
    
    Delta.neg.t = [Delta.neg.t; per(pksneg(:))];
    Delta.pos.t = [Delta.pos.t; per(pkspos(:))];
    Delta.neg.map = [Delta.neg.map; -seg(pksneg,:)];
    Delta.pos.map = [Delta.pos.map; seg(pkspos,:)];
    end
 
    
end
thr = 1*std(teigall);
delpos = LocalMinima(-teigall, MinISI, -thr)*eSampleRate/ResampleTo;
delneg = LocalMinima(teigall, MinISI, -thr)*eSampleRate/ResampleTo;

%lab = [repmat('neg',length(Delta.neg.t),1); repmat('pos',length(Delta.pos.t),1)];
lab = [ones(length(Delta.neg.t),1); ones(length(Delta.pos.t),1)*2; ones(length(delneg),1)*3; ones(length(delpos),1)*4];
res = [Delta.neg.t(:); Delta.pos.t(:); delneg(:); delpos(:)];

[res, ind ] = sort(res);
lab = lab(ind);
MakeEvtFile(res,[FileBase '.del.evt'],lab,eSampleRate);

