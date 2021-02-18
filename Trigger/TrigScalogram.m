%function [TrigSpecAv,TrigSpecStd,TrigWave,f,t,s] = TrigScalogram(X,T,Window,Fs, FreqRange)
%window in seconds
function [TrigSpecAv,TrigSpecStd,TrigWave,f,t,s] = TrigScalogram(X,T,varargin)

[Window,Fs, FreqRange] = DefaultArgs(varargin,{2,1250,[1 50]});

% nSamples = 2*Window*Fs+1;
Per = [-Window*Fs:Window*Fs];
if (min(size(X))>1) 
    error('cant do more than one dimension!');
end
X=X(:); %deal with one dimension so far
maxTime = size(X,1);
nColumns = size(X,2);
MyTriggers = T(find(T>Window*Fs & T<maxTime-Window*Fs));
nTriggers = length(MyTriggers);

TrigSpecAv =[];
TrigSpecStd =[];
TrigWave = [];

if nTriggers==0
    f=[];t=[];s=[];
    return
end

h =waitbar(0,'Wait ..');
for i=1:nTriggers
    [yo,f,t,s]=Wavelet(X(MyTriggers(i)+Per),FreqRange,Fs,6);
    pow = abs(yo).^2;
    if (i==1)
        TrigSpecAv = pow;
        TrigSpecStd = pow.^2;
        TrigWave = yo;
    else
        TrigSpecAv = TrigSpecAv+pow;
        TrigSpecStd = TrigSpecStd+pow.^2;
        TrigWave = TrigWave+yo;
    end
    waitbar(i/nTriggers,h);
    
end
close(h);
TrigSpecAv = TrigSpecAv/nTriggers;
TrigSpecStd = TrigSpecStd/nTriggers - TrigSpecAv.^2;
TrigWave = TrigWave/nTriggers;
%keyboard
% yo=reshape(yo,[length(f),nTriggers,nFFT]);

% AllSpec1=[];
% for i=1:nColumns
%     AllSpec1(:,:,i) = repmat(squeeze(AllSpec(:,i,i)),[1, length(t)]);
% end
%AllSpec = repmat(AllSpec1,[1, length(t) ,nColumns]);
%keyboard
% spec = log(TrigSpec) - log(AllSpec1);
if (nargout ==0)
   % figure
   for i=1:nColumns
       subplot(nColumns,1,i);
        imagelog(t,f,log(squeeze(TrigSpecAv(:,:,i)))');
       % set(gca,'YLim',[0 200]);
    end
end
