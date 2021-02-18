function out = PhaseSynchronization(FileBase,varargin)
%function PhaseSynchronization(FileBase,fMode, ElLoc, FreqRange,State,FrThreshold,Window)
% the name is wrong - it computes, in fact, unit - lfp coherence (and other measures) for each
% channel for each cell.
[fMode, ElLoc, FreqRange,States,FrThreshold,Window] = DefaultArgs(varargin, ...
    {'compute','c',[30 150], {'REM','RUN','SWS'},3,2^8});

Par = LoadPar([FileBase '.xml']);
MinPeriod = 5; %seconds for theta periods selection
if strcmp(ElLoc,'all')
    El = [1:Par.nElecGps];
elseif isstr(ElLoc(1))
    El = find(strcmp(Par.ElecLoc,ElLoc));
else
    El = ElLoc(1);
end
if ischar(States); States = {States}; end;
nSt = length(States);

[RepCh IfTetrode Info] = RepresentChan(FileBase);
if Info.nChInShank==8 IfTetrode =1; end
nRepCh = length(RepCh);
eFs = 1250;

switch fMode
    case 'compute'
        [AllRes,AllClu,Map] = LoadCluRes(FileBase,El);
        if isempty(AllRes)
            fprintf('no units in %s in this file \n',ElLoc);
            return;
        end
        AllRes = round(AllRes/16)+1;
        if ~isstr(ElLoc)
            myi = find(Map(:,2)==ElLoc(1) & Map(:,3)==ElLoc(2));
            AllRes  = AllRes(AllClu==myi);
            AllClu = ones(length(AllRes),1);
        end
        sc=1;
        %loop through States
        for s=1:nSt
            %Period = SelectStates(FileBase,States{s},MinPeriod);
            if ~FileExists([FileBase '.sts.' States{s}])
                fprintf('No state %s\n',States{s});
                continue;
            end

            Period = load([FileBase '.sts.' States{s}]);
            PerLen = round(sum(diff(Period,1,2)/1250));
            fprintf('Processing state %s of duration %d seconds\n',States{s},PerLen);

            [Res ind]=SelPerDiscr(AllRes,Period,1,1);
            Clu=AllClu(ind);

           
            fr = FiringRate(Res,Clu,[],1250);
            gfr = find(fr>FrThreshold);

            coh=[]; ph=[];powf=[];powu=[];
            % h = waitbar(0,'Please wait...');

            for ii=1:nRepCh
                eeg = LoadBinary([FileBase '.eeg'],RepCh(ii),Par.nChannels,Period);
                % waitbar(ii/nRepCh,h);
                for jj=1:length(gfr);
                    myRes = Res(Clu==gfr(jj));
                    if length(myRes)==0 continue; end
                    [y, out(sc).f,  phi, yerr, phierr, phloc,pow] = mtptchd(eeg,myRes,...
                        ones(length(myRes),1),2^(nextpow2(Window)+1),1250,2^nextpow2(Window),[],3,[],[],FreqRange,[],[],[],2);
                    out(sc).powf(:,ii) = sq(y(:,1,1));
                    out(sc).powu(:,jj) = sq(y(:,2,2));
                   % out(sc).powu(:,ii,jj) = sq(pow(:,1,2,:));
                    out(sc).coh(:,ii,jj) = sq(y(:,1,2));
                    out(sc).ph(:,ii,jj) = sq(phi(:,1,2));
                    out(sc).plph(:,ii,jj) = sq(angle(phloc(:,1,2)));
                    out(sc).ploc(:,ii,jj) = sq(abs(phloc(:,1,2)));
                    out(sc).yerr(:,ii,jj) = sq(yerr(:,1,2,1));
                end
            end
            %close(h);
            out(sc).ElClu = Map(gfr,:);
            out(sc).FirRate = fr(gfr);
            out(sc).State = States{s};
            sc= sc+1;
        end %loop through states
        if nargout<1
            save([FileBase '.' mfilename '.mat'],'out');
        end

    case 'fix'
        load([FileBase '.' mfilename '.mat']);
        for i=length(out):-1:1
         if isempty(out(i).f) & isempty(out(i).powu)
             out(i)=[];
         end
        end
        save([FileBase '.' mfilename '.mat'],'out');
%             out(1).f=out(2).f; 
%         else
%             return
%         end
        if isfield(out,'UseChannels') 
            out = rmfield(out,{'UseChannels','UseChannelsi','sig'});
            save([FileBase '.' mfilename '.mat'],'out');
        else
            return
        end
  %      keyboard
    case 'check'
       load([FileBase '.' mfilename '.mat']);
        if isfield(out,'sig')
            out=1;
        else
            out=0;
        end
        
        
    case {'display', 'select'}
                load([FileBase '.' mfilename '.mat']);
          %  global OutArgs;
        %  Map = MapSilicon(x,Channels,ChInShank)
        Channels = RepCh;
        Map = SiliconMap(Par);
        ChMap = Map.GridCoord(Channels,:);
        CluLoc = ClusterLocation(FileBase,El);
        for s=1:length(out)
            fprintf('%s\n',out(s).State);
            % GammaComponents(FileBase,'display',out(s).State);
            signum=1;
            figure(3232);clf
            for jj=1:size(out(s).coh,3)
                mycell = find(CluLoc(:,1)==out(s).ElClu(jj,1) & CluLoc(:,2)==out(s).ElClu(jj,2));
                if isempty(mycell)
                    continue;
                end
                if IfTetrode
                    for ll=1:length(Par.AnatGrps)
                        mychan1 = CluLoc(mycell,3);
                        if ismember(mychan1,Par.AnatGrps(ll).Channels+1)
                            mychan = RepCh(ll);
                            mychani = ll;
                        end
                    end
                else
                    mychan = CluLoc(mycell,3);
                    [dummy mychani] = min(abs(RepCh-mychan));
                end

                for k=1:Info.nShanks
                    if sum(Info.AllChannels{k}==mychan)
                        ych = find(Info.AllChannels{k}==mychan);
                        xch = k;
                        break;
                    end
                end
                if sum(out(s).powu(:,jj))==0
                    fprintf('did not fire enough\n');
                    continue;
                end
                subplot(5,2,[1 3])
                imcoh = sq(out(s).coh(:,:,jj));
                %                    imcoh(:,mychani) = NaN;
                imcohbutch =  imcoh(:,setdiff([1:size(imcoh,2)],mychani));
                imagesc(out(s).f,[1:nRepCh],imcoh');hold on
                caxis([min(imcohbutch(:)) max(imcohbutch(:))]);
                colorbar('NorthOutside');
                title('coherences'); ylabel('channels');
                hl = Lines([],mychani,'w',[],2);
                hold off

                subplot(5,2,5);
                %plot(out(s).f,sq(out(s).powu(:,jj)));axis tight
                hpy = plotyy(out(s).f,sq(out(s).powu(:,jj)),out(s).f,log10(sq(out(s).powf(:,mychani))));
                set(hpy(1),'XLim',[out(s).f(1),out(s).f(end)]);    set(hpy(2),'XLim',[out(s).f(1),out(s).f(end)]);
                title(['unit power - ' num2str(out(s).ElClu(jj,:)) ' ' ....
                    num2str(out(s).FirRate(jj)) ' Hz, State ' States{s}]);
                while 1
                    subplot(5,2,[1 3])
                    [x y b] = ginput(1);
                    if b==1
                        [dummy fi] = min(abs(out(s).f-x));
                        %[dummy chi] = min(abs([1:size(out(s).coh,3)]-y));
                        chi = round(y); chi = max(chi,1); chi = min(chi,size(out(s).coh,2));
                        
                        %get the location of clicked channel
                        for k=1:Info.nShanks
                            if sum(Info.AllChannels{k}==Channels(chi))
                                yclick = find(Info.AllChannels{k}==Channels(chi));
                                xclick = k;
                                break;
                            end
                        end
                        
                        mycoh = sq(out(s).coh(fi,:,jj));
                        subplot(3,2,2)
                        MapSilicon(mycoh,Par,[],0);
                        hold on
                        hpt = plot(xch, ych, 'ko','MarkerSize',20);
                        hpt1 = plot(xclick, yclick, 'kx','MarkerSize',20);
                        
                        hold off
                        if isfield(out(s), 'ploc')
                            myploc = sq(out(s).ploc(fi,:,jj));
                        else
                            myploc = 1;
                        end
                        subplot(3,2,4)
                        MapSilicon(myploc,Par,[],1);
                        hold on
                        hpt = plot(xch, ych, 'ko','MarkerSize',20);
                        hpt1 = plot(xclick, yclick, 'kx','MarkerSize',20);
                        hold off

                        mypow = sq(out(s).powf(fi,:));
                        myph = unwrap(sq(out(s).ph(fi,:,jj)));
                        subplot(3,2,6)
                        %MapSilicon(mypow,Par,[],1);
                        MapSilicon(myph,Par,16); 
                        hold on
                        hpt = plot(xch, ych, 'ko','MarkerSize',20);
                        hpt1 = plot(xclick, yclick, 'kx','MarkerSize',20);
                        hold off

                        subplot(5,2,7)
                        confplot(out(s).f,sq(out(s).coh(:,chi,jj)),sq(out(s).yerr(:,chi,jj)));
                        axis tight; ylim([0 max(ylim)]);
                        % plot(out(s).f,sq(out(s).coh(:,chi,jj)));axis
                        % tight
                        %hold on
                        % plot(out(s).f,sq(out(s).yerr(:,chi,jj)),'r--');
                        %hold off
                        title(['coherence to channel ' num2str(RepCh(chi))]);

                        subplot(5,2,9)
                        plot(out(s).f,unwrap(sq(out(s).ph(:,chi,jj)))*180/pi);axis tight
                        phax = ylim;
                        set(gca,'YTick',SelectPeriods([-180*5:90:180*5],phax,'d',1));
                        phaxlab = str2num(get(gca,'YTickLabel'));
                        phaxlab = mod(phaxlab,360);
                        set(gca,'YTickLabel',num2str(phaxlab));
                        title(['phase shift to channel ' num2str(RepCh(chi))]);
                        
                    elseif b==3
                        break;
                    elseif b==2
                        keyboard
                    else % presses the key  = get params of the coherence
                        if strcmp(fMode,'select')
                            out(s).sig(signum).CellInd = jj;
                            out(s).sig(signum).ElClu = out(s).ElClu(jj,:);
                            subplot(5,2,[1 3])
                            fprintf('select freq. ranges on the left plot\n');
                            [frx dummy] = ginput(2);
                            %get the left-right boundaries
                            [m fri] = min(abs(dotdot(frx,'-',out(s).f(:)')),[],2);
                            fri = sort(fri);
                            out(s).sig(signum).UseFreq = out(s).f(fri);
                            out(s).sig(signum).UseFreqi = [fri(1):fri(end)];

                            if ~IfTetrode
                                fprintf('Choose how to select channels (left - clicking,right -clustering)\n');
                                [xch ych b] = ginput(1);
                            else
                                b=1;
                                fprintf('use clicking - this is tetrode\n');
                            end
                            
                            switch b
                                case 1 %just click on good channel
                                    uch = zeros(length(Channels),1);
                                    pth = [];
                                    while 1
                                        subplot(3,2,2)
                                        fprintf('enter a channel : \n');
                                        [xch ych b] = ginput(1);
                                        [dd chi] = min(abs(ChMap(:,1)-xch).^2+abs(ChMap(:,2)-ych).^2);
                                        switch b
                                            case 1
                                                uch(chi)=1;
                                                fprintf('Channel %d is on\n', Channels(chi));
                                            case 3
                                                uch(chi)=0;
                                                fprintf('Channel %d is off\n', Channels(chi));
                                            case 2
                                                break;
                                            otherwise
                                                break;
                                        end
                                    end
                                    if IfTetrode
                                        out(s).sig(signum).UseChannels = Channels(uch==1);
                                        out(s).sig(signum).UseChannelsi = find(uch);
                                    else
                                        selCh = find(uch==1); selCh = sort(selCh);
                                        out(s).sig(signum).UseChannels =  Channels([selCh(1):selCh(end)]);
                                        out(s).sig(signum).UseChannelsi = [selCh(1):selCh(end)];
                                    end
                                case 3 % cluster channels on the map
                                    fprintf('now cluster the channels\n');
                                    in = ClusterPoints(ChMap,0);
                                    out(s).sig(signum).UseChannels = Channels(in);
                                    out(s).sig(signum).UseChannelsi = find(in);
                            end
                            signum=signum+1;
                            break;
                        else
                            prteps([mfilename '_' FileBase '_' num2str(out(s).ElClu(jj,:))]);
                            %return;
                        end
                    end
                end %while


            end
        end %loop over states
        if strcmp(fMode,'select')
            save([FileBase '.' mfilename '.mat'],'out');
        end
        
    case 'group'
        load([FileBase '.' mfilename '.mat']);
        in = out; clear out;
        Channels = RepCh;
        Map = SiliconMap(Par);
        ChMap = Map.GridCoord(Channels,:);
        CluLoc = ClusterLocation(FileBase,El);

        xyz = ChannelsXYZ(FileBase,RepCh,Par);
        if ~isfield(in,'sig') 
            out = [];
            return; 
        end
            
        cnt=1;
        for s=1:length(in)
            if isempty(in(s).sig) continue; end
            for k=1:length(in(s).sig)
                mycell = find(CluLoc(:,1)==in(s).sig(k).ElClu(1) & CluLoc(:,2)==in(s).sig(k).ElClu(2));
                if isempty(mycell)
                    continue;
                end
                if IfTetrode
                    mychan1 = CluLoc(mycell,3);
                    for ll=1:length(Par.AnatGrps)
                        if ismember(mychan1,Par.AnatGrps(ll).Channels+1)
                            mychan = RepCh(ll);
                            mychani = ll;
                        end
                    end
                else
                    mychan = CluLoc(mycell,3);
                    [dummy mychani] = min(abs(RepCh-mychan));
                end

                if IfTetrode & length(in(s).sig(k).UseChannelsi)==1
                    goodch = in(s).sig(k).UseChannelsi;
                elseif IfTetrode
                    goodch = setdiff(in(s).sig(k).UseChannelsi,mychani);
                else
                    goodch = setdiff(in(s).sig(k).UseChannelsi,mychani+[-1:1]); %remove 3 channels around 
                end
                goodf = in(s).sig(k).UseFreqi;
                f = in(s).f(goodf);
                jj = in(s).sig(k).CellInd;
                mycoh =  sq(in(s).coh(goodf,goodch,jj));
                myphi =  sq(in(s).ph(goodf,goodch,jj));
                if 0
                    %subplot(211)
                    figure(9991);clf
                    imagesc(f,[1:length(goodch)],mycoh'); colorbar
                    %                subplot(212)
                    figure(9992);clf
                    imagesc(f,[1:length(goodch)],myphi'*180/pi);
                    caxis([-180 180]);
                    CircColormap;
                    colorbar
                    waitforbuttonpress
                end
                out(cnt).State = in(s).State;
                out(cnt).FileBase = FileBase;
                out(cnt).El= in(s).sig(k).ElClu(1);
                out(cnt).Clu= in(s).sig(k).ElClu(2);
                out(cnt).Fr = in(s).FirRate(jj);
                mcoh = median(mycoh);
                [maxcoh maxchi] = max(mcoh);
                [out(cnt).MaxCoh maxfi] = max(mycoh(:,maxchi));
                out(cnt).Phi = myphi(maxfi,maxchi);
                out(cnt).Freq = f(maxfi);
                p = polyfit(f,unwrap(myphi(:,maxchi)),1);
                out(cnt).FreqSlope =p(1);
                out(cnt).TimeLag = p(1)/2/pi*1000;
                out(cnt).SigGap = out(cnt).MaxCoh - sq(in(s).yerr(goodf(maxfi),goodch(maxchi),jj));
                out(cnt).PhLoc = sq(in(s).ploc(goodf(maxfi),goodch(maxchi),jj));
                
                %cellxy = [CluLocXY.ML(mycell) CluLocXY.AP(mycell)];
                %maxxy = 
                MaxCohCh = RepCh(goodch(maxchi));
                CellCh = mychan;
                xyzMaxCoh = ChannelsXYZ(FileBase,MaxCohCh,Par);
                xyzCell = ChannelsXYZ(FileBase,CellCh,Par);
                                
                out(cnt).DistXY = sqrt(sum((xyzMaxCoh(1:2)-xyzCell(1:2)).^2));
                out(cnt).DistXYZ = sqrt(sum((xyzMaxCoh-xyzCell).^2));
                out(cnt).DistZ =abs(xyzMaxCoh(3)-xyzCell(3));
                cnt=cnt+1;
            end   
                

        end % end of state loop
        out = CatStruct(out);
end





