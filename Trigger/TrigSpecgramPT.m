%function [spec, f, t] = TrigSpecgramPT(X,res, clu, T,Window,nWindows,Fs)
%  X - matrix of continuous signal,
%  res and clu - spike times and cluster indexes,
%  T  - trigger times ( same sampling as the X and res)
%  Window - size of spectral windows(time bins)
%  nWidnows - number of spectral windows : nWindows*Window will be the size of the spectrogeram. 
function [spec, f,t] = TrigSpecgramPT(X,res, clu, T,Window,nWindows,Fs)
 
if (nargin<3 |isempty(Window))
    Window = 2^10;
end

if nargin<4 | isempty(nWindows) nWindows = 10; end
if length(Window)==1
    nSamples = 2*Window*nWindows;
    WindowSize = Window;
    Window = [Window*nWindows Window*nWindows];
else
    %don't know what this does :))
    nSamples = sum(Window);
    Window  = Window*nWindows;
end

if (nargin<5 | isempty(Fs))
    Fs = 1250;
end

%X=X(:); %deal with one dimension so far
maxTime = size(X, 1);
nColumns = size(X,2);
MyTriggers = T(find(T>Window(1) & T<maxTime-Window(2)));
BlockSize = floor(2000000/nSamples); % memory saving parameter
nTriggers = length(MyTriggers);

if nTriggers==0
    spec =[];f=[];t=[];
    return
end


% XTriggered = zeros(nTriggers,nSamples);
% go through triggers in groups of BlockSize to save memory
% for Block = 1:ceil(nTriggers/BlockSize)
% 	BlockTriggers = MyTriggers(1+(Block-1)*BlockSize:min(Block*BlockSize,nTriggers));
% 	nBlockTriggers = length(BlockTriggers);
% 	
%     TimeOffsets = repmat(-Window(1)+1:Window(2), nBlockTriggers, 1);
% 	TimeCenters = repmat(T(BlockTriggers), 1, nSamples);
% 	TimeIndex = TimeOffsets + TimeCenters;
%     Waves = X(TimeIndex);
% 	%Waves = reshape(Waves, [nBlockTriggers, nSamples, nColumns]);
%     XTriggered(1+(Block-1)*BlockSize:min(Block*BlockSize,nTriggers),:) = Waves;
% end

%WrapXTrig = reshape(XTriggered',[1,nTriggers*nSamples]);
nFFT = nSamples/nWindows;

%[AllSpec, f] = mtcsd(X(1:min(Fs*300,maxTime),:),nFFT,Fs,WindowSize,0);

TrigSpec =[];%zeros(length(f),nWindows*2-2,nColumns);
for i=1:nTriggers
    [yo, f,t]=mtptchd(X(MyTriggers(i)+[(-Window(1)+1):Window(2)],:),res, clu, nFFT,Fs,WindowSize,WindowSize/2); %,NW,Detrend,nTapers);
    if (i==1)
        TrigSpec =yo;
    else
        TrigSpec = TrigSpec + yo;
    end
end
%keyboard
% yo=reshape(yo,[length(f),nTriggers,nFFT]);
TrigSpec = TrigSpec/nTriggers;
AllSpec1=[];
for i=1:nColumns
    AllSpec1(:,:,i) = repmat(squeeze(AllSpec(:,i,i)),[1, length(t)]);
end
%AllSpec = repmat(AllSpec1,[1, length(t) ,nColumns]);
%keyboard
spec = log(TrigSpec) - log(AllSpec1);
if (nargout ==0)
   % figure
   for i=1:nColumns
       subplot(nColumns,1,i);
        imagesc(t,f,squeeze(spec(:,:,i)));
        set(gca,'YLim',[0 200]);
    end
end
