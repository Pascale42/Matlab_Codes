% function OutArgs = Gamma2LFP(FileBase, fMode, SubSet, SpkGps, CluLocOnly, State)
% 
% FileBase : Root name of the experiment in a string
% fMode: 'compute' or 'display'
% SubSet : Name it as you wish or 'All' ;
% SpkGps : enter all the SpikeGroup numbers belonging to the Subset if chosen, ex for AnatGp2: [6:10] or [6 7 8 9 10]
% CluLocOnly : run the function only on the clusters recorded by this list of channels (from 1)
% State : 'THE, 'SWS'....
% 

function OutArgs = Gamma2LFP(FileBase, fMode, Name, LfpCh, Win, State)

switch fMode
    case 'compute'
        Par = LoadPar([FileBase '.xml']);


        load([FileBase '.DetectGammaBurst.mPFC' State '.mat'], 'GBtime');
        gam = GBtime; % for a given state

           %%% Load EEG
        
        EEG=LoadBinary([FileBase '.eeg'],  LfpCh, Par.nChannels);
        SWS = load([FileBase '.sts.SWS']);
        EEG=SelectPeriods(EEG(:,:), SWS, 'c', 1);
        
        % Triger Average
        
        win=round(Par.lfpSampleRate*Win); % fenetre pour le triggered average
        % Faire l'axe x
        Srange = [-win(1):win(end)];
        Trange = linspace(-win(1),win(end),length(Srange));
        Trange=Trange/Par.lfpSampleRate;
        
        
        % retirer les indices de bords
        a=find(rips-win >0,1,'first');
        b=find(rips +win < size(EEG,1), 1, 'last');
        rips=rips(a:b);
        clear a b
 
     
        %%%  Triggered average Ripples LFP
        
        lfp=NaN(length(rips), ((win*2)+1));
        for n=1:length(rips)   % for each ripple
            lfp(n,:)=EEG([rips(n)-win : rips(n)+win]) ;
        end
        Av=mean(lfp,1);
        clear lfp
        
        OutArgs.Trange=Trange;
        OutArgs.Av=Av;
        
        save([FileBase '.' mfilename '.' Name  '.mat'], 'OutArgs')
        
    case 'display'
        
        load([FileBase '.' mfilename '.' Name  '.mat'], 'OutArgs')
        figure('name',[mfilename ' - ' Name],'NumberTitle','off')
        plot(OutArgs.Trange,OutArgs.Av); axis tight; hold on
        line([0 0],[min(OutArgs.Av) max(OutArgs.Av)])
        title(['Triggered LFP ' Name ' with Ripples']);
end
