% function Ripples2LFP(FileBase, fMode, Name, LfpCh, Win)
%
% fMode :: 'compute' or 'display' or 'dispRec'
% Name :: as you wish
% LfpCh :: Channel number for the LFP Trigger
% Win :: window for triggering, in seconds

function Ripples2LFP(FileBase, fMode, Name, varargin)

Par = LoadPar([FileBase '.xml']);

[LfpCh, Win] = DefaultArgs(varargin,{1, 0.25});

switch fMode
    case 'compute'
        
        %%% Loading Ripples
        
        if ~FileExists([FileBase '.spw.mat'])
            fprintf('no spws');
            return;
        end
        load([FileBase '.spw.mat']);
        rips = Rips.t; clear Rips
        
        
        %%% Load EEG
        
        EEG=LoadBinary([FileBase '.eeg'], LfpCh,  Par.nChannels);
%         SWS = load([FileBase '.sts.SWS']);
%         EEG=SelectPeriods(EEG(:,:), SWS, 'c', 1);
        
        %%% Triger Average
        
        win=round(Par.lfpSampleRate*Win); % fenetre pour le triggered average
        % Faire l'axe x
        Srange = -win(1):win(end);
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
            lfp(n,:)=EEG(rips(n)-win : rips(n)+win) ;
        end
        Av=mean(lfp,1);
        
        OutArgs.lfp=lfp;
        clear lfp
        OutArgs.Trange=Trange;
        OutArgs.Av=Av;
        
        save([FileBase '.' mfilename '.' Name  '.mat'], 'OutArgs')
        
    case 'display'
        
        load([FileBase '.' mfilename '.' Name  '.mat'], 'OutArgs')
        figure('name',[mfilename ' - ' Name],'NumberTitle','off')
        plot(OutArgs.Trange,OutArgs.Av); axis tight; hold on
        line([0 0],[min(OutArgs.Av) max(OutArgs.Av)], 'Color', [1 0 0])
        title(['Triggered LFP ' Name ' with Ripples']);
        
    case 'dispRec'
        
        load([FileBase '.' mfilename '.' Name  '.mat'], 'OutArgs')
        figure('name',[mfilename ' - ' Name],'NumberTitle','off')
        title(['Triggered LFP ' Name ' with Ripples']);
        for n=1: length(OutArgs.lfp)
            plot(OutArgs.Trange,OutArgs.lfp(n,:)); axis tight; hold on
            line([0 0],[min(OutArgs.lfp(n,:)) max(OutArgs.lfp(n,:))], 'Color', [1 0 0])
            waitforbuttonpress; clf
        end
        
end

