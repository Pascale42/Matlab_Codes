%
% a function to compute and display the modulation of cell firing by oscill phase
%
% function: OutArgs = PhaseModulation(FileBase, varargin :: State, fMode, Channel, SubSet,SpkGps, FreqRange, CluLocOnly, Overwrite)
%
% varargin are by default: 'THE','compute', 15, 'All', [], [3 6], [], 1
%
% fMode: 'compute'; 'display';  'displaynosort'; 'dispcluloc'; 'displaySingle'
% State: THE, SWS, ...
% FreqRange: [3 6] Hz for THE; [0.1 2] for SWS; [25 55] for GAM;
% Freely moving: [6 10] for THE; [0.1 5] for SWS
% Channel: takes eeg from channel 1; (starts at 0 in the xml!!!)
% SubSet : Name it as you wish or 'All' ;
% SpkGps : enter all the SpikeGroup numbers belonging to the Subset if chosen, ex for AnatGp2: [6:10] or [6 7 8 9 10]
% CluLocOnly : run the function only on the clusters recorded by this list of channels (from 1)
% ResCoef: downsample again the data, 1 by default
% Overwrite: 0=no  1=yes
% SaveFig: 0=no  1=yes
%
% The output is a structure OutArgs containing info on phase and units, as
% well as Map, with the correspondance between units # and shanks #



function OutArgs = PhaseModulation(FileBase, varargin)

[State, fMode, Channel, SubSet,SpkGps, FreqRange, CluLocOnly, Overwrite] = ...
    DefaultArgs(varargin,{'THE','compute', 15, 'All',[],[3 6], [], 1});


switch fMode
    case 'compute'
        
        
        %%% Basic checkings and loading FileBase.CluRes.mat
        
        if exist([FileBase '.PhaseModulation' State '.mat'], 'file')>0 && ~Overwrite
            disp('Already computed!');return
        end
        
        %%% Loading units

        
        if strcmp(SubSet, 'All');
            if ~FileExists([FileBase '.CluRes.mat'])
                [T,G,Map, Par]=LoadCluRes(FileBase);
            else
                load([FileBase '.CluRes.mat']);
            end
        else
            [T,G,Map,Par]=LoadCluRes(FileBase, SpkGps,0, CluLocOnly);
        end
        
        %%% Loading the states
        
        if strmatch(State, 'Gamma') %#ok<*MATCH2>
            Stat= input('Which State?', 's');
            Structure = input('Gamma in which Structure?', 's');
            load([FileBase '.DetectGammaBursts.' Structure '.' Stat '.mat'])
            STA=GammaBursts;
        else
            if ~exist([FileBase '.sts.' State], 'file')
                disp('No sts file for this state');return
            end
            STA = load([FileBase '.sts.' State]);
        end
        
        T = round(T/(Par.SampleRate/Par.lfpSampleRate));  % spikes are sampled @ 20000 or 32556 Hz whereas eeg is 1250Hz
        [T, it]=unique(T); % remove possible double detections;  12/20/2016
        G=G(it); clear it
        
        
        %%% Loading eeg and getting good segments and the spikes in the segments
        
        eeg  = LoadBinary([FileBase '.eeg'], Channel, Par.nChannels);
        [eeg, orind]= SelectPeriods(eeg(:),STA,'c',1);
                
        [Tsta, Indsta] = SelPerDiscr2(T, STA);
        Gsta = G(Indsta);
        
        [uClu, ~, NewClu] = unique(Gsta);
        % OutArgs.Map = [[1:length(uClu)]' Map(uClu,2:3)]; removed 03/2015
        x= ismember(Map(:,1), uClu);
        OutArgs.Map = Map(x,:); clear x
        
                
        %%% Filtering the eeg in the Frequency range and getting continuous phase
        
        feeg = ButFilter(eeg, 2, FreqRange/(Par.lfpSampleRate/2),'bandpass');
        hilb = hilbert(feeg);
        ph = angle(hilb);
        
        
        %%% Rayleigh statistic for phase modulation
        
        phstats = PhaseOfTrain(ph(ismember(orind,Tsta)),NewClu,[],[],1,0,max(unique(NewClu))); % orind implemented 11/24/2016
        for n=1:length(uClu) % in case of an empty cluster
            if isempty(phstats(:,n).Clu);
                phstats(:,n).Clu=n;
                phstats(:,n).pval=0;
                phstats(:,n).R=0;
                phstats(:,n).th0=0;
                phstats(:,n).phconfint=[0;0];
                phstats(:,n).mu=0;
                phstats(:,n).k=0;
                phstats(:,n).phhist=zeros(1,73);
                phstats(:,n).phbin=[-pi:pi/18:3*pi];
            end
        end
        
        
        %%% Saving the output as OutArgs structure array
        
        OutArgs.cphstats  = CatStruct(phstats);
        OutArgs.ChannelUsed = Channel;
        
        if strcmp(SubSet,'All')
            save([FileBase '.' mfilename State '.mat'],'OutArgs');
        else
            save([FileBase '.' mfilename State '.' SubSet '.mat'],'OutArgs');
        end
        
        
        
        
        
    case 'display'   % Classical display mode
        
        
        if strcmp(SubSet,'All')
            load([FileBase '.' mfilename State '.mat']);
        else
            load([FileBase '.' mfilename State '.' SubSet '.mat']);
        end
        
        
        figure('name',['Phase Modulation - ' SubSet],'NumberTitle','off')
        
        if length(OutArgs.cphstats.Clu) == 1
            subplot(121)
            [xl, yl]=pol2cart(1, 1); c = compass(xl,yl); set(c, 'Visible', 'off'); hold on % R lim 1
            polar(OutArgs.cphstats.th0, OutArgs.cphstats.R,'o');
            title(['Phase of ' num2str(State) ' on which cells fire']);
            
            subplot(122)
            bar(OutArgs.cphstats.phbin*180/pi, OutArgs.cphstats.phhist)
            xlim([-200 600]); xlabel(['THE phase, deg']);
            title('Phase modulation'); axis tight
        else
            
            subplot(131)
            [xl, yl]=pol2cart(1, 1); c = compass(xl,yl); set(c, 'Visible', 'off'); hold on % R lim 1
            polar(OutArgs.cphstats.th0, OutArgs.cphstats.R,'o');
            title(['Phase of ' num2str(State) ' on which cells fire']);
            
            subplot(132)
            [ s, si ] = sort(OutArgs.cphstats.th0);
            imagesc(OutArgs.cphstats.phbin(:,1)*180/pi, [1:size(OutArgs.cphstats.phhist,2)],...
                unity(OutArgs.cphstats.phhist(:,si))');
            xlabel([num2str(State) ' phase, deg']); ylabel('Cell index: Shank #, Cluster #');
            title([ num2str(State) ' Phase modulation']);
            set(gca,'YTick',[1:size(OutArgs.Map,1)]); colormap('jet')
            set(gca,'YTickLabel',num2str(OutArgs.Map(si,2:3)));
            
            subplot(133)
            rose(OutArgs.cphstats.th0,60); title('Mean phases distribution');
            
            
            figure('name',['Phase Modulation - ' SubSet ' - individual neurons'],'NumberTitle','off')
            
            nclu=length(OutArgs.cphstats.Clu);
            for n=1:nclu
                subplotfit(n, nclu);  bar(OutArgs.cphstats.phbin(:,1)*180/pi, OutArgs.cphstats.phhist(:,n))
                xlim([-200 600]); xlabel(['phase, deg']); axis tight
                title(['Neuron # ' num2str(OutArgs.Map(n,(2:3)))]);
            end
         
        end
        
  

    case 'displaySignif'   % Display distinction between significantly entrained and not neurons 
        
           
        if strcmp(SubSet,'All')
            load([FileBase '.' mfilename State '.mat']);
        else
            load([FileBase '.' mfilename State '.' SubSet '.mat']);
        end
        
        x=find(OutArgs.cphstats.pval <= 0.5); xx=find(OutArgs.cphstats.pval > 0.5);
        
        figure('name',['Phase Modulation Signif Neurons - ' SubSet],'NumberTitle','off') 
            subplot(131)
            polar(OutArgs.cphstats.th0(x), OutArgs.cphstats.R(x),'om'); hold on
            polar(OutArgs.cphstats.th0(xx), OutArgs.cphstats.R(xx),'ob');
            title(['Phase of ' num2str(State) ' on which cells fire']);
           
            
            subplot(132)
            [ s, si ] = sort(OutArgs.cphstats.th0(x));
            imagesc(OutArgs.cphstats.phbin(:,1)*180/pi, [1:size(OutArgs.cphstats.phhist(:,x),2)],...
                unity(OutArgs.cphstats.phhist(:,si))');
            xlabel([num2str(State) ' phase, deg']); ylabel('Cell index: Shank #, Cluster #');
            title([ num2str(State) ' Phase modulation']);
            set(gca,'YTick',[1:size(OutArgs.Map(x,:),1)]); colormap('jet')
            set(gca,'YTickLabel',num2str(OutArgs.Map(si,2:3)));
            
            subplot(133)
            rose(OutArgs.cphstats.th0(x),60); title('Mean phases distribution');
            
            
            figure('name',['Phase Modulation - ' SubSet ' - individual Signif neurons'],'NumberTitle','off')
            nclu=length(x);
            for n=1:nclu
                subplotfit(n, nclu);  bar(OutArgs.cphstats.phbin(:,1)*180/pi, OutArgs.cphstats.phhist(:,x(n)))
                xlim([-200 600]); xlabel(['phase, deg']); axis tight
                title(['Neuron # ' num2str(OutArgs.Map(x(n),(2:3)))]);
            end
  
        
        
    case 'displaySingle'   % Display a single neuron
        
        if strcmp(SubSet,'All')
            load([FileBase '.' mfilename State '.mat']);
        else
            load([FileBase '.' mfilename State '.' SubSet '.mat']);
        end
        
        clu=input('Enter the Cluster ID you want to plot (as [shank clu]) = ');
        [~, iclu]=Intersection(OutArgs.Map(:,[2 3]), clu);
        
        if isempty(iclu)
            disp(['The neuron ' num2str(clu) ' does not fire during ' State ' state'])
            return
        end
        
        figure('name',['Phase Modulation - ' SubSet],'NumberTitle','off')
        subplot(121)
        [xl, yl]=pol2cart(1, 1); c = compass(xl,yl); set(c, 'Visible', 'off'); hold on % R lim 1
        polar(OutArgs.cphstats.th0(iclu), OutArgs.cphstats.R(iclu),'o');
        title(['Phase of ' num2str(State) ' on which neuron ' num2str(clu) ' fires']);
        
        subplot(122)
        bar(OutArgs.cphstats.phbin(:,1)*180/pi, OutArgs.cphstats.phhist(:, iclu))
        xlim([-200 600]); xlabel(['THE phase, deg']);
        title('Phase modulation'); axis tight
        
        
    case 'displaynosort'   % Display no sorting in subplot 122
        
        
        if strcmp(SubSet,'All')
            load([FileBase '.' mfilename State '.mat']);
        else
            load([FileBase '.' mfilename State '.' SubSet '.mat']);
        end
        
        figure('name',['Phase Modulation - ' SubSet],'NumberTitle','off')
        
        subplot(121)
        [xl, yl]=pol2cart(1, 1); c = compass(xl,yl); set(c, 'Visible', 'off'); hold on % R lim 1
        polar(OutArgs.cphstats.th0, OutArgs.cphstats.R,'o');
        title(['Phase of ' num2str(State) ' on which cells fire']);
        
        subplot(122)
        imagesc(OutArgs.cphstats.phbin(:,1)*180/pi, [1:size(OutArgs.cphstats.phhist,2)],...
            unity(OutArgs.cphstats.phhist)');
        xlabel([num2str(State) ' phase, deg']); ylabel('Cell index: Shank #, Cluster #');
        title([ num2str(State) ' Phase modulation']);
        set(gca,'YTick',[1:size(OutArgs.Map,1)]);
        set(gca,'YTickLabel',num2str(OutArgs.Map(:,2:3)));
        
        
        figure('name',['Phase Modulation - ' SubSet ' - individual neurons'],'NumberTitle','off')
        
        nclu=length(OutArgs.cphstats.Clu);
        for n=1:nclu
            subplotfit(n, nclu);  bar(OutArgs.cphstats.phbin(:,1)*180/pi, OutArgs.cphstats.phhist(:,n))
            xlim([-200 600]); xlabel(['THE phase, deg']); axis tight
            title(['Phase modulation Neuron # ' num2str(OutArgs.Map(n,(2:3)))]);
        end
        
        
    case 'dispcluloc'   % Display with respect to ClusterLocation
        
        Par=LoadPar([FileBase '.xml']);
        if strcmp(SubSet,'All')
            load([FileBase '.' mfilename State '.mat']);
        else
            load([FileBase '.' mfilename State '.' SubSet '.mat']);
        end
        
        if FileExists([FileBase '.cluloc'])
            cluloc=  load([FileBase '.cluloc']);
        else
            cluloc=ClusterLocation(FileBase,[],0);
        end
        
        CluLoc=[];
        
        for n=1:length(OutArgs.Map)
            for j = 1:length(cluloc)
                if OutArgs.Map(n,2) == cluloc(j,1) && OutArgs.Map(n,3) == cluloc(j,2)
                    CluLoc=[CluLoc;cluloc(j,:)];
                end
            end
        end
        
        if  ismember(SpkGps, [1:5])
            x=[Par.AnatGrps(1).Channels+1]';
        end
        if ismember(SpkGps, [6:10])
            x=[Par.AnatGrps(2).Channels+1]';
        end
        
        
        CluLoc=[CluLoc zeros(size(CluLoc,1),1)];
        for n=1:size(CluLoc,1)
            CluLoc(n,4)=find(x == CluLoc(n,3));
        end
        [ s, si ] = sort(CluLoc(:,4));
        
        
        figure('name',['Phase Modulation - ' SubSet],'NumberTitle','off')
        
        subplot(121)
        [xl, yl]=pol2cart(1, 1); c = compass(xl,yl); set(c, 'Visible', 'off'); hold on % R lim 1
        polar(OutArgs.cphstats.th0, OutArgs.cphstats.R,'o');
        title(['Phase of ' num2str(State) ' on which cells fire']);
        
        subplot(122)
        imagesc(OutArgs.cphstats.phbin(:,1)*180/pi, [1:size(OutArgs.cphstats.phhist,2)],...
            unity(OutArgs.cphstats.phhist(:,si))');
        xlabel([num2str(State) ' phase, deg']); ylabel('Cell index: Shank #, Cluster #');
        title([ num2str(State) ' Phase modulation']); colorbar
        set(gca,'YTick',[1:size(OutArgs.Map,1)]);
        set(gca,'YTickLabel',num2str(CluLoc(si,1:3)));
        
         
        figure('name',['Phase Modulation - ' SubSet ' - individual neurons'],'NumberTitle','off')
         
         nclu=length(OutArgs.cphstats.Clu);
         for n=1:nclu
         subplotfit(n, nclu);  bar(OutArgs.cphstats.phbin(:,1)*180/pi, OutArgs.cphstats.phhist(:,n))
            xlim([-200 600]); xlabel(['THE phase, deg']); axis tight
            title(['Phase modulation Neuron # ' num2str(OutArgs.Map(n,(2:3)))]);
         end
         
end


