


function LFP = ConcatGamma(FileBase, WhereGam, Channels, Name, State, Win)




Par = LoadPar([FileBase '.xml']);
load([FileBase '.DetectGammaBursts.' WhereGam '.' State '.mat'], 'GBtime');
gam = GBtime; % for a given state


%%% Load EEG

EEG=LoadBinary([FileBase '.eeg'], Channels, Par.nChannels);
STA = load([FileBase '.sts.' State]);
EEG=SelectPeriods(EEG(:,:), STA, 'c', 1);



%%% Fenetre

win=round(Par.lfpSampleRate*Win);
% Srange = [-win(1):win(end)];
% Trange = linspace(-win(1),win(end),length(Srange));
% Trange=Trange/Par.lfpSampleRate;


% retirer les indices de bords
a=find(gam-win >0,1, 'first');
b=find(gam +win < size(EEG,1), 1, 'last');
gam=gam(a:b);
clear a b


%%%  Triggered Gamma LFP

LFP=[];
    for n=1:length(gam)   
        LFP=[LFP; EEG([gam(n)-win : gam(n)+win],:)];
    end




save([FileBase '.' mfilename '.' Name '.' State '.mat'], 'LFP')






