function out = ThGamMod(FileBase,varargin)
% out = ThGamMod(FileBase,fMode, States,LoadAll, ThBins,GamBins,Shuffle,PhaseChannels)
% LoadAll = 1 if you have enough memory to load all eeg channels into
% mememory (could be long for sws/run)
% ThBins and GamBins - nF x 2 matrix, beginning and end of freq. ranges to
% filter the eeg to compute phase or power (via hilbert transform)
% if ThBins = 0 then will use phase from thpar.mat
%ThBinsDef =dotdot([2.5:0.5:13]','+',[-1 1]);
ThBinsDef = 0;
%GamBinsDef=dotdot([30:10:40]','+',[-5 5]);
GamBinsDef=dotdot([30:3:200]','+',[-5 5]);
%ShuffleDef = struct('Type','shift','nShuffle',200,'MaxShift',125); % shift by 1 sec (after 10 times resampling);
ShuffleDef = 0;
%PhChannelDef = load([FileBase '.eegseg.par'])+1;
PhChannelDef  = '.usechan';
[fMode, States,LoadAll, ThBins,GamBins, Shuffle, PhChannel ] = ...
    DefaultArgs(varargin,{'compute',{'REM','RUN'},0, ThBinsDef, GamBinsDef, ShuffleDef,PhChannelDef});

if isstr(PhChannel) & FileExists([FileBase PhChannel])
    PhChannel = load([FileBase PhChannel])
end
%profile on;

% ThBins=ThBins(4:5,:);
% GamBins=GamBins(10:11,:);

Par = LoadPar([FileBase '.xml']);
[Channels,IfTetrode, Info] = RepresentChan(Par,'silicon');
%Channels =Channels(1); %%%% DEBUGGING!
PhChannelInd = find(ismember(Channels,PhChannel));
nStates = length(States);
MinPeriod = 1;

switch fMode
    case 'compute'
        wcnt=1;
        for w = 1:length(States) %loop through states
            if FileExists([FileBase '.sts.' States{w}])
                AllPeriod = load([FileBase '.sts.' States{w}]);
                Period = AllPeriod(find(diff(AllPeriod,1,2)>Par.lfpSampleRate*MinPeriod),:);
            else
                fprintf('No %s periods. Return empty output\n',States{w});
                continue;
            end

            if isempty(Period);
                fprintf('No %s periods. Return empty output\n',States{w});
                continue;
            end

            fprintf('Computing : State %s, \n',States{w});
            out(wcnt).State = States{w};
            out(wcnt).Channels = Channels;
            out(wcnt).fPh = mean(ThBins,2);
            out(wcnt).fPow = mean(GamBins,2);
            
            if ThBins==0
                load([FileBase '.thpar.mat'],'ThPh');
                ThPh = SelectPeriods(ThPh,Period,'c',1);
                if LoadAll==1
                    eeg = LoadBinary([FileBase '.lfp'],Channels, Par.nChannels, Period);
                    out(wcnt).thgam = PowerPhasePairs(ThPh, 0, eeg,  GamBins/625,10,@PowerModulation,Shuffle);
                else
                    loadChannel = @(ChInd) LoadBinary([FileBase '.lfp'],Channels(ChInd), Par.nChannels, Period);
                    eeg = struct('FunHandle',loadChannel,'nPow',length(Channels));
                    out(wcnt).thgam = PowerPhasePairs(ThPh,0,eeg ,  GamBins/625,10,@PowerModulation,Shuffle);
                end
            else
                if LoadAll==1
                    eeg = LoadBinary([FileBase '.eeg'],Channels, Par.nChannels, Period);
                    out(wcnt).thgam = PowerPhasePairs(eeg(:,PhChannelInd), ThBins/625, eeg,  GamBins/625,10,@PowerModulation,Shuffle);
                else
                    loadChannel = @(ChInd) LoadBinary([FileBase '.eeg'],Channels(ChInd), Par.nChannels, Period);
                    eeg = struct('FunHandle',loadChannel,'nPow',length(Channels));
                    out(wcnt).thgam = PowerPhasePairs(eeg.FunHandle(PhChannelInd), ThBins/625,eeg ,  GamBins/625,10,@PowerModulation,Shuffle);
                end
            end
            wcnt=wcnt+1;
        end

        save([FileBase '.' mfilename '.mat'],'out');

    case 'display'
        load([FileBase '.' mfilename '.mat']);
        if ndims(out(1).thgam.Ramp)==2 error('not yet working for one phase channel frequency only');end
       % keyboard
        myFig = figure(667);
        %menuDisp = uimenu(myFig, 'Label','States','HandleVisibility','off');

        nStates = length(out);
        States ={};
        for ii=1:nStates
            States{ii} = out(ii).State;
       %     funH(ii) = uimenu(menuDisp,'Label',out(ii).State, 'Callback',@figupdate);
        end
%        keyboard
        %set initial values
        nChannels = length(out(1).Channels);
        nfPow = length(out(1).fPow);
        nfPh = length(out(1).fPh);
        Map = SiliconMap(Par);
        ChMap = Map.GridCoord(Channels,:);
        w=1;
        PP='Pow';
        [~, fPowi] = min(abs(out(1).fPow-50));
        [~, fPhi] = min(abs(out(1).fPh-8));
        Chi = 1;
        PhChi = 1;
        AmpType='Ramp'; %1 -Ramp, 2- MI
        SigType='Rzscore'
        getpos=0;

        set(myFig, 'WindowButtonDownFcn', @figupdate);
        set(myFig, 'KeyPressFcn', @varupdate);
        myh=[];
        figupdate();
%    keyboard
end
    function figupdate(varargin)
        figure(myFig);
        
        set(myFig,'NumberTitle','off','Name',['ThGamMod ' States{w}]);

        if getpos
            [axN figN axH figH] = WhereIClicked;
            axN = find(myh==axH);
            whatbutton = get(gcf,'SelectionType');
            mousecoord = get(gca,'CurrentPoint');
            x=mousecoord(1,1); y = mousecoord(1,2);
            switch axN
                case 1 
                    [~, fPhi] = min(abs(out(1).fPh-x));
                    [~, fPowi] = min(abs(out(1).fPow-y));
                case 2
                    [~, Chi] = min(abs([1:nChannels]-round(x)));
                    [~, fPowi] = min(abs(out(1).fPow-y));
                case 3
                    [~, fPhi] = min(abs(out(1).fPh-x));
                    [~, Chi] = min(abs([1:nChannels]-round(y)));
                case 4
                    [dd fPowi] = min(abs(out(1).fPow-y));
                case 5
                    [~, fPhi] = min(abs(out(1).fPh-y));
                case {6,7}  % Map for all channels at given freqs
                    [dd Chi] = min((ChMap(:,1)-x).^2+(ChMap(:,2)-y).^2);
                    fprintf('Current channel %d\n',Channels(Chi)); 
                case 8
                    [dd fPowi] = min(abs(out(1).fPow-x));
                case 9
                    [dd fPhi] = min(abs(out(1).fPh-x));
                    
            end
        end
        
        %now do theplots
        nrow=3; ncol=3;
       
        clf
        myh(1) = subplot2(nrow,ncol,1,1); % 1: amplitude fPow x fPh
        if strcmp(PP,'Ph')
            imagescnan({out(w).fPh,out(w).fPow, sq(out(w).thgam.Rth(:,:,Chi,PhChi))'*180/pi}, [-180 180], 1,1);
            axis xy;  title('Theta Phase');
        else
            imagesc(out(w).fPh,out(w).fPow, sq(out(w).thgam.(AmpType)(:,:,Chi,PhChi))'); axis xy; colorbar; title(AmpType);
        end
        xlabel('fPh'); ylabel('fPow');
        hold on
        plot(out(w).fPh(fPhi), out(w).fPow(fPowi),'ko');
        LineCross(out(w).fPh(fPhi), out(w).fPow(fPowi),'w',[],2);
       
        
        
        myh(2) = subplot2(nrow,ncol,1,2); % 2: amplitude fPow x ChInd
        if strcmp(PP,'Ph')
            imagescnan({[1:nChannels],out(w).fPow, sq(out(w).thgam.Rth(fPhi,:,:,PhChi))*180/pi}, [-180 180], 1,1);axis xy
            title('Theta Phase');
        else
            imagesc([1:nChannels],out(w).fPow, unity(sq(out(w).thgam.(AmpType)(fPhi,:,:,PhChi)))); axis xy; colorbar;  title(AmpType);
        end
        ylabel('fPow'); xlabel('Channels');
        hold on
        plot(Chi, out(w).fPow(fPowi),'ko');
        LineCross(Chi,out(w).fPow(fPowi),'w',[],2);
       

        myh(3) = subplot2(nrow,ncol,2,1); % 3: amplitude fPh x ChInd
        if strcmp(PP,'Ph')
            imagescnan({out(w).fPh,[1:nChannels], sq(out(w).thgam.Rth(:,fPowi,:,PhChi))'*180/pi}, [-180 180], 1,1);axis xy
            title('Theta Phase');
        else
            imagesc(out(w).fPh,[1:nChannels], (sq(out(w).thgam.(AmpType)(:,fPowi,:,PhChi)))'); axis xy; colorbar
        end
        xlabel('fPh'); ylabel('Channels');
        hold on
        plot(out(w).fPh(fPhi),Chi, 'ko');
        LineCross(out(w).fPh(fPhi),Chi,'w',[],2);
        title(AmpType);
                
        myh(4) =subplot2(nrow,ncol,1,3); % 4: phase density vs fPow
        phbin= out(1).thgam.phbins(:,1,1,1,1);
        imagesc([phbin(:); phbin(:)+2*pi]*180/pi,out(1).fPow,...
            [sq(out(w).thgam.pow_dens(:,fPhi,:,Chi,PhChi)); sq(out(w).thgam.pow_dens(:,fPhi,:,Chi,PhChi))]'); 
        axis xy
        ylabel('fPow'); xlabel('Theta Phase');
        hold on
        Lines([], out(w).fPow(fPowi), 'k');
        
        myh(5) =subplot2(nrow,ncol,2,3); % 5: phase density vs fPh
        phbin= out(1).thgam.phbins(:,1,1,1,1);
        imagesc([phbin(:); phbin(:)+2*pi]*180/pi,out(1).fPh,...
            [sq(out(w).thgam.pow_dens(:,:,fPowi,Chi,PhChi)); sq(out(w).thgam.pow_dens(:,:,fPowi,Chi,PhChi))]'); 
        axis xy
        ylabel('fPh'); xlabel('Theta Phase');
        hold on
        Lines([], out(w).fPh(fPhi), 'k');
        
        if 1
        myh(6) =subplot2(nrow,ncol,2,2); % 6: amplitude map
        MapSilicon(sq(out(w).thgam.(AmpType)(fPhi,fPowi,:,PhChi)), [], Par,[],0,[],0,[],1);
%        MapSilicon(-log10(eps+sq(out(w).thgam.Rpval(fPowi,:))), [], Par,[],0,[],0,[],1);
        xlabel('ML'); ylabel('DV');
        hold on
        plot(ChMap(Chi,1),ChMap(Chi,2),'ko');
        LineCross(ChMap(Chi,1),ChMap(Chi,2),'w',[],2);

        
        myh(7) =subplot2(nrow,ncol,3,2); % 7: phase map 
        %sigmask = double(sq(out(w).thgam.Rpval(fPowi,:))>0.01);
        %sigmask(sigmask==1)=NaN;
        %sigmask(sigmask==0)=1;
        sigmask = ones(nChannels,1);
        MapSilicon(sq(out(w).thgam.Rth(fPhi,fPowi,:,PhChi))*180/pi.*sigmask, [], Par,[],0,[-180 180],1,[],1);
        xlabel('ML'); ylabel('DV');
        hold on
        plot(ChMap(Chi,1),ChMap(Chi,2),'ko');
        LineCross(ChMap(Chi,1),ChMap(Chi,2),'w',[],2);
        title('phase');
        end
        myh(8) =subplot2(nrow,ncol,3,1); % 8: modulation by phase reference
        plot(out(w).fPow, sq(out(w).thgam.(AmpType)(fPhi,:,Chi,:)));
        legend(num2str([1:size(out(1).thgam.Rth,4)]'));
        hold on
        Lines(out(1).fPow(fPowi),[],'k');
        xlabel('fPow');axis tight
        title(['current phase from ' num2str(PhChi)]);
        
        myh(9) =subplot2(nrow,ncol,3,3); % 9: modulation by phase reference
        plot(out(w).fPh, sq(out(w).thgam.(AmpType)(:,fPowi,Chi,:)));
        legend(num2str([1:size(out(1).thgam.Rth,4)]'));
        hold on
        Lines(out(1).fPh(fPhi),[],'k');
        xlabel('fPh');
        axis tight
        
%          myh(10) =subplot2(nrow,ncol,4,3);
%          if strcmp(PP,'Ph')
%             imagescnan({out(w).fPh,[1:nChannels], sq(out(w).thgam.Rth(:,fPowi,:,PhChi))'*180/pi}, [-180 180], 1,1);axis xy
%             title('Theta Phase');
%         else
%             imagesc(out(w).fPh,[1:nChannels], sq(out(w).thgam.(AmpType)(:,fPowi,:,PhChi))'); axis xy; colorbar
%         end
%         
%         xlabel('fPh'); ylabel('Channels');
%         hold on
%         plot(out(w).fPh(fPhi),Chi, 'ko');
%         LineCross(out(w).fPh(fPhi),Chi,'w',[],2);
%         title(AmpType);
%       
%         myh(6) =subplot(nrow,ncol,ncol+1); % 6: sig map
%         MapSilicon(sq(out(w).thgam.(SigType)(fPowi,:)), [], Par,[],0,[],0,[],1);
%         xlabel('ML'); ylabel('DV');
%         hold on
%         plot(ChMap(Chi,1),ChMap(Chi,2),'ko');
%         LineCross(ChMap(Chi,1),ChMap(Chi,2),'w',[],2);
%         title(SigType)
%         
%         myh(7) =subplot(nrow,ncol,ncol*2+1); % 7: histogram at current channel, fPow and fPh
%         %phe = linspace(-pi,pi,size(out(w).thgam.pow_dens,1)+1); phbin = (phe(1:end-1)+phe(2:end))/2;
%         phbin= out(1).thgam.phbins(:,1,1,1,1);
%         barcirc(phbin, sq(out(w).thgam.pow_dens(:,fPowi,Chi)));
%         ylim(sum(ylim.*[0 1]).*[0.5 1.2]); 
%         xlabel('Phase');
%          
         
%         myh(9) = subplot(nrow,ncol,ncol*3); % 9: phase fPow x ChInd
%         imagescnan({[1:nChannels], out(w).fPow, out(w).thgam.Rth*180/pi}, [-180 180], 1,1); axis xy; 
%         ylabel('fPow'); xlabel('Channels');
%         hold on
%         plot(Chi, out(w).fPow(fPowi),'ko');
%         LineCross(Chi,out(w).fPow(fPowi),'w',[],2);
        
        
        if getpos==0; getpos=1; end
%        fstepPow=1; fstepPh=1;
    end

    function varupdate(varargin)
         whatkey = get(myFig,'CurrentCharacter');
         switch double(whatkey)
             case double('s')
                 % change state value
                 if nStates==1
                     fprintf('just one state');
                 else
                     if w==1 w=2; else w=1; end
                     fprintf('Current state %s\n',States{w});
                     getpos=0;
                     figupdate();
                 end
             case double('1')
                 PhChi = 1;
                  getpos=0;
                 figupdate();
             case double('2')
                 PhChi = 2;
                  getpos=0;
                 figupdate();
             case double('3')
                 PhChi = 3;
                  getpos=0;
                 figupdate();
              case double('4')
                 PhChi = 4;
                  getpos=0;
                 figupdate();
             case double('d')
                 %change display
                 if strcmp(AmpType,'Ramp') 
                     AmpType='MI';
                     SigType='MIzscore';
                 else
                     AmpType='Ramp';
                     SigType='Rzscore';
                 end
                 getpos=0;
                 figupdate();
                 
             case double('p')
                 if strcmp(PP,'Pow')
                     PP = 'Ph';
                 else
                    PP = 'Pow';
                    
                 end
               
                 figupdate();
                 
             case 29
                 fPowi = min(fPowi+1,length(out(1).fPow));figupdate();
             case 28
                 fPowi = max(1,fPowi-1);figupdate();
                 
             case double('p')
                 %print
                 prteps(['ThGamMod -' datestr(now)]);
                 
             case double('q');
                    set(myFig, 'WindowButtonDownFcn',[]);
                    set(myFig, 'KeyPressFcn', []);
                    return;
             
         end
    end
end
    
