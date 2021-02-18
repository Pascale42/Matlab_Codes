%
% a function to compute and display the modulation of cell firing by oscill phase
%
% function: OutArgs = PhaseModulation2(FileBase, varargin)
%
% varargin: State, fMode, Channel, SubSet, SpkGps, FreqRange, ResCoef, Overwrite, SaveFig
%  are by default: 'THE','compute', 15, 'All', [], [3 6], 1, 1, 0
%
% fMode: 'compute'; 'display';  'displaynosort'; 'dispcluloc';
% State: THE, SWS, GAMTHE, GAMSWS
% FreqRange: [3 6] Hz for THE; [0.1 2] for SWS; [25 55] for GAM;
% Channel: takes eeg from channel 1; (starts at 0 in the xml!!!)
% SubSet : 'AnatGp1', 'AnatGp2', 'All'  corresponds to the spikes of the first AnatGroup (Hippocampus), or of the second, or all of them
% SpkGps : enter all the SpikeGroup numbers belonging to the Subset if chosen, ex for AnatGp2: [6:10] or [6 7 8 9 10]
% ResCoef: downsample again the data, 1 by default
% Overwrite: 0=no  1=yes
% SaveFig: 0=no  1=yes
%
% The output is a structure OutArgs containing info on phase and units, as
% well as Map, with the correspondance between units # and shanks #



function OutArgs = PhaseModulation2(FileBase, varargin)

[State, fMode, Channel, SubSet,SpkGps, FreqRange, ResCoef, Overwrite, SaveFig] = ...
    DefaultArgs(varargin,{'THE','compute', 15, 'All',[],[3 6], 1, 1, 0});


switch fMode
    case 'compute'
        
        
        %%% Basic checkings and loading FileBase.CluRes.mat
        
        if exist([FileBase '.PhaseModulation' State '.mat'])>0 & ~Overwrite
            disp('Already computed!');return
        end
        
%         
%         if ~FileExists([FileBase '.CluRes.mat'])
%             Par = LoadPar([FileBase '.xml']);
%             [T,G,Map]=LoadCluRes(FileBase);
%         else
%             load([FileBase '.CluRes.mat']);
%         end
%         

if strcmp(SubSet, 'All');
    
    if ~FileExists([FileBase '.CluRes.mat'])
        [T,G,Map, Par]=LoadCluRes(FileBase);
    else
        load([FileBase '.CluRes.mat']);
    end
else
    
    [T,G,Map,Par]=LoadCluRes(FileBase, SpkGps,0);
    
end

        
        
        
        %%% Loading the states
        
        if strcmp(State, '9Hz');
            State='THE';
            
            if ~exist([FileBase '.sts.' State])
                disp('No sts file for this state');return
            end
            STA = load([FileBase '.sts.' State]);
            State='9Hz';
        else
            if ~exist([FileBase '.sts.' State])
                disp('No sts file for this state');return
            end
            STA = load([FileBase '.sts.' State]);
        end
        
        T = round(T/(Par.SampleRate/Par.lfpSampleRate));  % !!!! spikes are @ 32556 Hz whereas eeg is 1250Hz !!!!!
        
        %%% Loading eeg and getting good segments and the spikes in the segments
        
        eeg  = LoadBinary([FileBase '.eeg'], Channel, Par.nChannels);
        [eeg orind]= SelectPeriods(eeg(:),STA,'c',1);
        
        
        [Tsta Indsta] = SelPerDiscr(T, STA, 1, 1);
        Gsta = G(Indsta);
        
        [uClu dummy NewClu] = unique(Gsta); 
        clear dummy 
%         OutArgs.Map = [[1:length(uClu)]' Map(uClu,2:3)]; Changed with Andi
        OutArgs.Map = [uClu Map(uClu,2:3)];
        
        %%% Resampling the eeg: let's be lighter for slow rhythms
        
        eeg = resample(eeg,1,ResCoef);
        
        eSampleRate = Par.lfpSampleRate/ResCoef;
        
        Tsta = round(Tsta/ResCoef);
        Tsta(Tsta==0) = 1;
        Tsta(Tsta>length(eeg)) = length(eeg);
        
        
%         %%% Getting a SubSet of neurons or not
%         
%         switch SubSet
%             
%             case 'AnatGp1'
%                 
%                 clu=find(OutArgs.Map(:,2) <= SpkGps(end), 1,'last');               
% 
%                 clu=OutArgs.Map(clu);
%                 indx=find(Gsta<=clu);
%                 Gsta=Gsta(indx);
%                 Tsta=Tsta(indx);
%                 [uClu dummy NewClu] = unique(Gsta);
%                 clear dummy
%                 OutArgs.Map = [uClu Map(uClu,2:3)];
%                 
%                 
%             case 'AnatGp2'
%                 
%                 clu=find(OutArgs.Map(:,2) >= SpkGps(1), 1,'first');
%                 clu=OutArgs.Map(clu); % get the real cluster number
%                 
%                 indx=find(Gsta >= clu);
%                 Gsta=Gsta(indx);
%                 Tsta=Tsta(indx);
%                 [uClu dummy NewClu] = unique(Gsta);
%                 clear dummy
%                  OutArgs.Map = [uClu Map(uClu,2:3)];
%                 
%         end
        
        
        %%% Filtering the eeg in the Frequency range and getting continuous phase

        feeg = ButFilter(eeg, 2, FreqRange/(eSampleRate/2),'bandpass');
        hilb = hilbert(feeg);
        ph = angle(hilb);
        
        
        %%% Rayleigh statistic for phase modulation
        
        phstats = PhaseOfTrain(ph(Tsta),NewClu,[],[],1,0,max(unique(NewClu)));
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
        
        
        
        
        
    case 'display'
        
        
        if strcmp(SubSet,'All')
            load([FileBase '.' mfilename State '.mat']);
        else
            load([FileBase '.' mfilename State '.' SubSet '.mat']);
        end
        
        
        figure('name',['Phase Modulation - ' SubSet],'NumberTitle','off')
        
        subplot(131)
        [xl yl]=pol2cart(1, 1); c = compass(xl,yl); set(c, 'Visible', 'off'); hold on % R lim 1
        polar(OutArgs.cphstats.th0, OutArgs.cphstats.R,'o');
        title(['Phase of ' num2str(State) ' on which cells fire']);
        
        subplot(132)
        [ s si ] = sort(OutArgs.cphstats.th0);
        imagesc(OutArgs.cphstats.phbin(:,1)*180/pi, [1:size(OutArgs.cphstats.phhist,2)],...
            unity(OutArgs.cphstats.phhist(:,si))');
        xlabel([num2str(State) ' phase, deg']); ylabel('Cell index: Shank #, Cluster #');
        title([ num2str(State) ' Phase modulation']);
        set(gca,'YTick',[1:size(OutArgs.Map,1)]);
        set(gca,'YTickLabel',num2str(OutArgs.Map(si,2:3)));
        
        subplot(133)
        bar(OutArgs.cphstats.phbin(:,1)*180/pi, sum(OutArgs.cphstats.phhist,2))
        xlim([-200 600]); xlabel(['THE phase, deg']);
        title('THE Phase modulation');
        
        if SaveFig > 0;
            reportfig(gcf,['PhaseModulation_' State], 0, ['File ' FileBase ', State ' State],150);
        else
            
        end
        
        
        
        
        
    case 'displaynosort'
        
        
        if strcmp(SubSet,'All')
            load([FileBase '.' mfilename State '.mat']);
        else
            load([FileBase '.' mfilename State '.' SubSet '.mat']);
        end
        
        figure('name',['Phase Modulation - ' SubSet],'NumberTitle','off')
        
        subplot(131)
        [xl yl]=pol2cart(1, 1); c = compass(xl,yl); set(c, 'Visible', 'off'); hold on % R lim 1
        polar(OutArgs.cphstats.th0, OutArgs.cphstats.R,'o');
        title(['Phase of ' num2str(State) ' on which cells fire']);
        
        subplot(132)
        imagesc(OutArgs.cphstats.phbin(:,1)*180/pi, [1:size(OutArgs.cphstats.phhist,2)],...
            unity(OutArgs.cphstats.phhist)');
        xlabel([num2str(State) ' phase, deg']); ylabel('Cell index: Shank #, Cluster #');
        title([ num2str(State) ' Phase modulation']);
        set(gca,'YTick',[1:size(OutArgs.Map,1)]);
        set(gca,'YTickLabel',num2str(OutArgs.Map(:,2:3)));
        
        subplot(133)
        bar(OutArgs.cphstats.phbin(:,1)*180/pi, sum(OutArgs.cphstats.phhist,2))
        xlim([-200 600]); xlabel(['THE phase, deg']);
        title('THE Phase modulation');
        
        if SaveFig > 0;
            reportfig(gcf,['PhaseModulation' State 'nosort'], 0, ['File ' FileBase ', State ' State],150);
        end
        
        
        
        
        
    case 'dispcluloc'
        
        
        if strcmp(SubSet,'All')
            load([FileBase '.' mfilename State '.mat']);
        else
            load([FileBase '.' mfilename State '.' SubSet '.mat']);
        end
        
        
        cluloc=  load([FileBase '.cluloc']);
        CluLoc=[];
        
        for n=1:length(OutArgs.Map)
            for j = 1:length(cluloc)
                if OutArgs.Map(n,2) == cluloc(j,1) & OutArgs.Map(n,3) == cluloc(j,2)
                    CluLoc=[CluLoc;cluloc(j,:)];
                end
            end
        end
        
        
        x=[Par.AnatGrps(2).Channels+1]';
        CluLoc=[CluLoc zeros(size(CluLoc,1),1)];
        for n=1:size(CluLoc,1)
            CluLoc(n,4)=find(x(:,1) == CluLoc(n,3));
        end
        [ s si ] = sort(CluLoc(:,4));
        
        figure(30000)
        
        imagesc(OutArgs.cphstats.phbin(:,1)*180/pi, [1:size(OutArgs.cphstats.phhist,2)],...
            unity(OutArgs.cphstats.phhist(:,si))');
        xlabel([num2str(State) ' phase, deg']); ylabel('Cell index: Shank #, Cluster #');
        title([ num2str(State) ' Phase modulation']);
        set(gca,'YTick',[1:size(OutArgs.Map,1)]);
        set(gca,'YTickLabel',num2str(CluLoc(si,1:3)));
end













%-----------------------------------------------------------------------------------------

%ph = ProcessInChunks(feeg,{'hilbert','angle'},2^16,2^10); % let's chunk
%amp = abs(hilb);
% plot(unity([eeg(1:1250) feeg(1:1250) ph(1:1250)]))
% phn = MakeUniformDistr(ph);



%         subplot(132)
%         imagesc(OutArgs.phbin(:,1), [1:size(OutArgs.phhist,2)], unity(OutArgs.phhist)');
%         [ s si ] = sort(OutArgs.cphstats.th0);
%         imagesc(OutArgs.cphstats.phbin(:,1)*180/pi, [1:size(OutArgs.cphstats.phhist,2)],...
%             unity(OutArgs.cphstats.phhist)');
%         xlabel([num2str(State) ' phase, deg']); ylabel('Cell index: Shank #, Cluster #');
%         title([ num2str(State) ' Phase modulation']);
%         set(gca,'YTick',[1:size(OutArgs.Map,1)]);
%         set(gca,'YTickLabel',num2str(OutArgs.Map(:,2:3)));
%         subplot(133)
%         [ s si ] = sort(OutArgs.cphstats.th0);
%         %normalisation of the histograms
%         a=OutArgs.cphstats.phhist;
%         sa = sum(a);
%         repsa=repmat(sa,size(a,1),1);
%         norma = a./repsa;
%         imagesc(OutArgs.cphstats.phbin(:,1)*180/pi, [1:size(norma,2)],...
%             (norma(:,si))');
%         xlabel([num2str(State) ' phase, deg']); ylabel('Cell index: Shank #, Cluster #');
%         title([ num2str(State) ' Phase modulation, normalized']);
%         set(gca,'YTick',[1:size(OutArgs.Map,1)]);
%         set(gca,'YTickLabel',num2str(OutArgs.Map(si,2:3)));
