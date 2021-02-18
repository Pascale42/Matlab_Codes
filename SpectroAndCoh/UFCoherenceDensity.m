%
% Unit Field Multitaper Coherence Density
% function OutArgs = UFCoherenceDensity(FileBase, State, fMode, Channels, FreqRange, SubSet, SpkGps, CluLocOnly);
%
% SubSet : Name it as you wish or 'All' ;
% SpkGps : enter all the SpikeGroup numbers belonging to the Subset if chosen, ex for AnatGp2: [6:10] or [6 7 8 9 10]
% CluLocOnly : run the function only on the clusters recorded by this list of channels (from 1)
% State : 'THE, 'SWS'....
% SaveFig: 1 to save report
% fMode: 'compute' and 'display' modes
% Channels, ResCoef, FreqRange: Don't give if you want same channels as the
% ones computed in LoadFilEeg
% |-------------------------|
% |          |  unit- unit    |
% |          |       coh       |
% |_____|__________|
% | 3       |    unit-lfp     |
% | 2       |     coh         |
% | 1 2 3 |__________|
% | lfp            units      |





function OutArgs = UFCoherenceDensity(FileBase, State, fMode, Channels, FreqRange, SubSet, SpkGps, CluLocOnly, State2)

if nargin < 9
    State2 = '';
end



switch fMode
    
    case 'compute'
        
        
        Par = LoadPar([FileBase '.xml']);
        WinLengthSec=1;
        % loading State
        
        if strmatch(State, 'Gamma')
            load([FileBase '.DetectGammaBursts.mPFC.'  State2 '.mat'])
            STA=GammaBursts;
        else
            if ~exist([FileBase '.sts.' State], 'file')
                disp('No sts file for this state');return
            end
            STA = load([FileBase '.sts.' State]);
        end
        
        
        
        % loading all channels
        
        eeg = LoadBinary([FileBase '.eeg'], Channels, Par.nChannels);
        eeg = SelectPeriods(eeg(:,:),STA,'c',1);
        
        % loading Units;
        
        if strcmp(SubSet, 'All');
            if ~FileExists([FileBase '.CluRes.mat'])
                [T,G,Map, Par]=LoadCluRes(FileBase);
            else
                load([FileBase '.CluRes.mat']);
            end
        else
            [T,G,Map,Par]=LoadCluRes(FileBase, SpkGps,0, CluLocOnly);
        end
        T =  max(1,round(T*Par.lfpSampleRate/Par.SampleRate));
        [T, Ind] = SelPerDiscr(T, STA, 1, 1);
        G = G(Ind);
        [uClu, ~, NewClu] = unique(G);
%         OutArgs.Map = [uClu Map(uClu,2:3)];
        x= ismember(uClu, Map(:,1));
        OutArgs.Map = Map(x,:); clear x
        
        
        WinLengthSample = 2^round(log2(WinLengthSec*Par.lfpSampleRate));
        nFFT = 2*WinLengthSample;
        OutArgs.FreqRange = FreqRange;
        OutArgs.Channels= Channels;
                
        tic
        [OutArgs.ufchd.y, OutArgs.ufchd.f, OutArgs.ufchd.phi, OutArgs.ufchd.yerr, OutArgs.ufchd.phierr, OutArgs.ufchd.phloc, OutArgs.ufchd.pow]=...
            mtptchd(eeg, T, G, nFFT,Par.lfpSampleRate,WinLengthSample,[],3,'linear',[],FreqRange);
        toc
        save([FileBase '.' mfilename '.' SubSet '.' State '.' State2 '.mat'], 'OutArgs');
        
        
    case 'display'
        load([FileBase '.' mfilename '.' SubSet '.' State '.' State2 '.mat']);
        
        nChannels = length(OutArgs.Channels);
        
        figure('name',[ mfilename ' - ' SubSet ' - ' State],'NumberTitle','off')
        for n=1:nChannels
            subplotfit(n, nChannels);
            imagesc(OutArgs.ufchd.f,[1:length(OutArgs.Map)],(squeeze(OutArgs.ufchd.y(:,n,(1+nChannels):end))'));
            set(gca,'YTick',[1:size(OutArgs.Map,1)]); set(gca,'YTickLabel',num2str(OutArgs.Map(:,2:3)));
            xlabel('Hz'); title(['Units-Field (ch ' num2str(OutArgs.Channels(n)) ') coherence density for ' num2str(State)]);
            ylabel('Neurons (spkgroup, cluster#)')
            
            colorbar; colormap('Jet')
        end
        
        
    case 'disp2'
        load([FileBase '.' mfilename '.' SubSet '.' State '.' State2 '.mat']);
        
        nChannels = length(OutArgs.Channels);
        nClu = size(OutArgs.Map,1);
        
        figure('name',[ mfilename ' - ' SubSet ' vs mPFC lpf - ' State],'NumberTitle','off')
        for n=1:nClu
            subplotfit(n,nClu);
            plot(OutArgs.ufchd.f, OutArgs.ufchd.y(:,1,n+nChannels))
            axis tight
            title(num2str(OutArgs.Map(n,2:3)));
        end
        
        figure('name',[ mfilename ' - ' SubSet ' vs Hpc lpf - ' State],'NumberTitle','off')
        for n=1:nClu
            subplotfit(n,nClu);
            plot(OutArgs.ufchd.f, OutArgs.ufchd.y(:,2,n+nChannels))
            axis tight
                title(num2str(OutArgs.Map(n,2:3)));
        end
       
    case 'sortdisp'
        
        
load([FileBase '.' mfilename '.' SubSet '.' State  '.' State2 '.mat']);
load([FileBase '.G2U-UFC.' SubSet '.' State2 '.mat'],'SortNeu');
%                     load([FileBase '.' mfilename '.' SubSet '.' State  '.mat']);
%           load([FileBase '.G2U-UFC.' SubSet '.' State '.mat'],'SortNeu');
           
          OutArgs.Map=[OutArgs.Map SortNeu];
         nChannels = length(OutArgs.Channels);
        x1=find(OutArgs.Map(:,4) == 1);
        x2=find(OutArgs.Map(:,4) == 2);
        x3=find(OutArgs.Map(:,4) == 3);
        
%         UFC=squeeze(OutArgs.ufchd.y(:,1,3:end));
%         
        figure
        % LFP channel 1
        subplot 321
        imagesc(OutArgs.ufchd.f,[1:length(OutArgs.Map(x1,:))],squeeze(OutArgs.ufchd.y(:,1,x1+nChannels))');
        axis tight
        set(gca,'YTick',[1:size(OutArgs.Map(x1,:),1)]);
        set(gca,'YTickLabel',num2str(OutArgs.Map(x1,2:3)));
        xlabel('Hz'); colorbar; colormap('Jet');clim([0 1]);
        
        subplot 323
        imagesc(OutArgs.ufchd.f,[1:length(OutArgs.Map(x2,:))],squeeze(OutArgs.ufchd.y(:,1,x2+nChannels))');
        axis tight
        set(gca,'YTick',[1:size(OutArgs.Map(x2,:),1)]);
        set(gca,'YTickLabel',num2str(OutArgs.Map(x2,2:3)));
        xlabel('Hz'); colorbar; colormap('Jet');clim([0 1]);
        
        subplot 325
        imagesc(OutArgs.ufchd.f,[1:length(OutArgs.Map(x3,:))],squeeze(OutArgs.ufchd.y(:,1,x3+nChannels))');
        axis tight
        set(gca,'YTick',[1:size(OutArgs.Map(x3,:),1)]);
        set(gca,'YTickLabel',num2str(OutArgs.Map(x3,2:3)));
        xlabel('Hz'); colorbar; colormap('Jet');clim([0 1]);
        
        % LFP channel 2
        subplot 322
        imagesc(OutArgs.ufchd.f,[1:length(OutArgs.Map(x1,:))],squeeze(OutArgs.ufchd.y(:,2,x1+nChannels))');
        axis tight
        set(gca,'YTick',[1:size(OutArgs.Map(x1,:),1)]);
        set(gca,'YTickLabel',num2str(OutArgs.Map(x1,2:3)));
        xlabel('Hz'); colorbar; colormap('Jet');clim([0 1]);
        
        subplot 324
        imagesc(OutArgs.ufchd.f,[1:length(OutArgs.Map(x2,:))],squeeze(OutArgs.ufchd.y(:,2,x2+nChannels))');
        axis tight
        set(gca,'YTick',[1:size(OutArgs.Map(x2,:),1)]);
        set(gca,'YTickLabel',num2str(OutArgs.Map(x2,2:3)));
        xlabel('Hz'); colorbar; colormap('Jet');clim([0 1]);
        
        subplot 326
        imagesc(OutArgs.ufchd.f,[1:length(OutArgs.Map(x3,:))],squeeze(OutArgs.ufchd.y(:,2,x3+nChannels))');
        axis tight
        set(gca,'YTick',[1:size(OutArgs.Map(x3,:),1)]);
        set(gca,'YTickLabel',num2str(OutArgs.Map(x3,2:3)));
        xlabel('Hz'); colorbar; colormap('Jet');clim([0 1]);
        
        if ~exist('OutArgs.MaxCoh', 'var')
        OutArgs.MaxCoh=NaN(size(OutArgs.Map,1),2);
        OutArgs.MaxCoh(:,1)=squeeze(max(OutArgs.ufchd.y(:,1,(1+nChannels):end)));
        OutArgs.MaxCoh(:,2)=squeeze(max(OutArgs.ufchd.y(:,2,(1+nChannels):end)));
        save([FileBase '.' mfilename '.' SubSet '.' State '.' State2 '.mat'], 'OutArgs');
        end
        
end