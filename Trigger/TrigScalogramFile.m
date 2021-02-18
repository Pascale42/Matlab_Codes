% function [Spec, f, t] = TrigScalogramFile(FileBase, Trigger, Channels, Window,FreqRange, ResampleCoef,FileExt)
% Window, SpecWindow and SpecStep are in msec
% Trigger in sampling rate of the file
function [Spec, f, t] = TrigScalogramFile(FileBase, Trigger,varargin)

Par =LoadPar(FileBase);
[Channels, Window, FreqRange, ResampleCoef, FileExt] = DefaultArgs(varargin,...
{[1:16],  1000,   [1 100], 1, 'eeg'});
eFs =1250;
rFs = (eFs/ResampleCoef);
nSamples = round(eFs*Window/1000);
nChannels = length(Channels);
nTriggers = length(Trigger);
StartPoints = Trigger-round(nSamples/2);
segs = LoadSegs([FileBase '.' FileExt], StartPoints,nSamples,Par.nChannels,Channels,eFs,ResampleCoef,0); 
[f,s,p] = WaveletSampling(rFs,FreqRange);
Spec=zeros(size(segs,1),length(f),nChannels);
hw = waitbar(0,'computing wavlet scalogram ..');
for i=1:nTriggers
    for j=1:nChannels
        [wave, t, f]=Wavelet(sq(segs(:,j,i)),FreqRange,rFs); 
        Spec(:,:,j) = Spec(:,:,j) +abs(wave.^2);
        waitbar(i/nTriggers,hw);
    end
end
close(hw);
    

return
SpecWindowSmpls = [SpecWin SpecStep]/1000;
spec = [];
for ii=1:nChannels
    mysegs = sq(segs(:,ii,:));
    sz = size(mysegs);
    if Whiten>0
        if ii==1
            [mysegs AR] = WhitenSignal(reshape(mysegs,[],1),[],1);
            mysegs = reshape(mysegs,sz);
        else
            mysegs = WhitenSignal(reshape(mysegs,[],1),[],1, AR);
            mysegs = reshape(mysegs,sz);
        end
    end
    params = struct('tapers',Tapers,'pad',1,'Fs',Fs/ResampleCoef,'fpass',FreqRange,'err',0,'trialave',Acc);
    if Acc==1
        [spec(:,:,ii),t,f]=mtspecgramc(mysegs,SpecWindowSmpls,params);%Tapers,1,Fs/ResampleCoef,FreqRange,0,1);
    else
        [spec(:,:,:,ii),t,f]=mtspecgramc(mysegs,SpecWindowSmpls,params);%Tapers,1,Fs/ResampleCoef,FreqRange,0,0);
    end
end
t =t - nSamples/Fs/2;
if nargout<1
    figure
    % matrix is t x f x ch - let's reshape it to make display stack of
    % specgrams along freq. for all channels with common time
    sz = size(spec);
    Mat = reshape(spec,[sz(1) sz(2)*sz(3)]);
    nf = sz(2)*sz(3);
    
    imagesc(t,[1:nf],Mat);

end

