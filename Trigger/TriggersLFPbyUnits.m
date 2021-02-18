% 
% function TriggersLFPbyUnits(FileBase, fMode, Name, varargin :: LfpCh, State, Win, SpkGps, CluLocOnly)
% 
%
% fMode :: 'compute', 'display'
% Name ::  Name it as you wish;
% LfpCh :: LFP channel (from 1) (default = [])
% State :: 'THE', 'SWS' ... (default = [])
% Win :: window for Trigger Average in seconds (default = 1])
% SpkGps : enter all the SpikeGroup numbers, ex for AnatGp2: [6:10] or [6 7 8 9 10] (default = [])
% CluLocOnly : run the function only on the clusters recorded by this list of channels (from 1) (default = [])
%
% 4 mai 2015

function TriggersLFPbyUnits(FileBase, fMode, Name, varargin)

Par = LoadPar([FileBase '.xml']);

[LfpCh, State, Win, SpkGps, CluLocOnly] = DefaultArgs(varargin,{1, [], 1, [],[]});

switch fMode
    case 'compute'
        

        
        %%% Loading units
        
        [T,G,Map,Par]=LoadCluRes(FileBase, SpkGps,0, CluLocOnly);
        T =  max(1,round(T*Par.lfpSampleRate/Par.SampleRate)); % get spike times to the LFP sample rate
        
        %%% Load EEG

        EEG=LoadBinary([FileBase '.eeg'],  LfpCh, Par.nChannels);
        
        %%% if State
        
        if ~isempty(State)
            if ~exist([FileBase '.sts.' State],'file')
                disp('No sts file for this state');return
            end
            STA = load([FileBase '.sts.' State]);
            
            [T, Ind] = SelPerDiscr(T, STA, 1, 1);
            G = G(Ind);
            
            EEG=SelectPeriods(EEG(:,:), STA, 'c', 1);
        end
        
        % Triger Average
        
        win=round(Par.lfpSampleRate*Win); % fenetre pour le triggered average
        % Faire l'axe x
        Srange = [-win(1):win(end)];
        Trange = linspace(-win(1),win(end),length(Srange));
        Trange=Trange/Par.lfpSampleRate;
        
        
        % retirer les indices de bords
        a=find(T-win >0,1,'first');
        b=find(T +win < size(EEG,1), 1, 'last');
        T=T(a:b);
        G=G(a:b);
        clear a b
        
        
        % --------  Triggered average RE LFP--------
        Av=NaN(size(Map,1),((win*2)+1));
        uclu=unique(G);
        
        for n=1:length(uclu)   % for each neuron
            indx=find(G == uclu(n));
            Tot=NaN(length(indx), ((win*2)+1));
            t=T(indx);
            for m=1:length(indx)
                Tot(m,:)=EEG([t(m)-win : t(m)+win]) ;
            end
            Av(n,:)=mean(Tot,1);
            clear t indx
        end
        clear Tot
        
        OutArgs.uclu=uclu;
        OutArgs.Trange=Trange;
        OutArgs.Av=Av;
        OutArgs.Map=Map;
        
        save([FileBase '.' mfilename '.' Name  '.' State '.mat'], 'OutArgs')
        
    case 'display'
        
        load([FileBase '.' mfilename '.' Name '.' State '.mat'], 'OutArgs')
        
        N=ceil(sqrt(length(OutArgs.uclu)));
        
        figure('name',[ mfilename ' - ' Name ' - ' State],'NumberTitle','off')
        for n=1:length(OutArgs.uclu)
            subplot(N,N,n)
            plot(OutArgs.Trange,OutArgs.Av(n,:)); axis tight; hold on
            line([0 0],[min(OutArgs.Av(n,:)) max(OutArgs.Av(n,:))])
            title(['Triggered LFP ' Name ' unit ' num2str(OutArgs.Map(n,2:3))]);
        end
        clear N
end

