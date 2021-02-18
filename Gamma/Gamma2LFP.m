% function OutArgs = Gamma2LFP(FileBase, fMode, WhereGam, WhereLfp, State, LfpCh, Win)
% 
% FileBase : Root name of the experiment in a string
% fMode: 'compute' or 'display'
% WhereGam : Where Gamma was detected (as used in DetectGammaBursts.m)
% WhereLfp : Name it as you wish 
% State : 'THE, 'SWS'....
% LfpCh :  Channel(s) for LFP
% Win : window in seconds
% 

function OutArgs = Gamma2LFP(FileBase, fMode, WhereGam, WhereLfp, State, LfpCh, Win)

switch fMode
    case 'compute'
        Par = LoadPar([FileBase '.xml']);


        load([FileBase '.DetectGammaBursts.' WhereGam '.' State '.mat'], 'GBtime');
        gam = GBtime; % for a given state

           %%% Load EEG
        
        EEG=LoadBinary([FileBase '.eeg'], LfpCh, Par.nChannels);
        STA = load([FileBase '.sts.'  State]);
        EEG=SelectPeriods(EEG(:,:), STA, 'c', 1);
        
        %%% Triger Average
        
        win=round(Par.lfpSampleRate*Win);
        Srange = [-win(1):win(end)];
        Trange = linspace(-win(1),win(end),length(Srange));
        Trange=Trange/Par.lfpSampleRate;
        
        
        % retirer les indices de bords
        a=find(gam-win >0,1,'first');
        b=find(gam +win < size(EEG,1), 1, 'last');
        gam=gam(a:b);
        clear a b
 
     
        %%%  Triggered average Gamma LFP
        
        % 1 or more channels in EEG
        if size(EEG,2) == 1
            lfp=NaN(length(gam), ((win*2)+1));
            for n=1:length(gam)  
                lfp(n,:)=EEG([gam(n)-win : gam(n)+win]) ;
            end
            Av=mean(lfp,1);
            clear lfp
        else
            for ch=1:size(EEG,2)
                lfp=NaN(length(gam), ((win*2)+1));
                for n=1:length(gam)  
                    lfp(n,:)=EEG([gam(n)-win : gam(n)+win],ch) ;
                end
                Av(ch,:)= mean(lfp,1); 
                clear lfp
            end
        end
        
        
        OutArgs.Trange=Trange;
        OutArgs.Av=Av;
        
        save([FileBase '.' mfilename '.' WhereGam '.' WhereLfp  '.' State '.mat'], 'OutArgs')
        
    case 'display'
        
        load([FileBase '.' mfilename '.' WhereGam '.' WhereLfp  '.' State '.mat'], 'OutArgs');
        
        figure('name',[mfilename ' - ' WhereGam '.' WhereLfp    ' - ' State ],'NumberTitle','off')
        nch =size(OutArgs.Av,1);
        for n=1:nch
            subplotfit(n,nch);
            plot(OutArgs.Trange,OutArgs.Av(n,:)); axis tight; hold on
            line([0 0],[min(OutArgs.Av(n,:)) max(OutArgs.Av(n,:))], 'Color', [1 0 0])
        end
  
end
