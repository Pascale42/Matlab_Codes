function out = PhaseSyncrhonization(FileBase,varargin)
%function PhaseSyncrhonization(FileBase,fMode, ElLoc, FreqRange,State,FrThreshold,Window)
[fMode, ElLoc, FreqRange,States,FrThreshold,Window] = DefaultArgs(varargin, ...
    {'compute','c',[30 150], {'REM','RUN'},1,2^8});

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

[RepCh ifTetr Info] = RepresentChan(FileBase);
nRepCh = length(RepCh);
eFs = 1250;

switch fMode
    case 'compute'
        [AllRes,AllClu,nclu,sph,elclu] = ReadEl4CCG(FileBase,El);
        if isempty(AllRes)
            fprintf('no units in %s in this file \n',ElLoc);
            return;
        end
        AllRes = round(AllRes/16)+1;
        if ~isstr(ElLoc)
            myi = find(elclu(:,1)==ElLoc(1) & elclu(:,2)==ElLoc(2));
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

            eeg = LoadBinary([FileBase '.eeg'],RepCh,Par.nChannels,Period);

            fr = FiringRate(Res,Clu,[],1250);
            gfr = find(fr>FrThreshold);

            coh=[]; ph=[];powf=[];powu=[];
           % h = waitbar(0,'Please wait...');

            for ii=1:nRepCh
               % waitbar(ii/nRepCh,h);
                for jj=1:length(gfr);
                    myRes = Res(Clu==gfr(jj));
                    if length(myRes)==0 continue; end
                    [y, out(s).f,  phi, yerr, phierr, phloc] = mtptchd(eeg(:,ii),myRes,...
                        ones(length(myRes),1),2^(nextpow2(Window)+1),1250,2^nextpow2(Window),[],2,[],[],FreqRange,[],[],[],2);
                    out(sc).powf(:,ii) = sq(y(:,1,1));
                    out(sc).powu(:,jj) = sq(y(:,2,2));
                    out(sc).coh(:,ii,jj) = sq(y(:,1,2));
                    out(sc).ph(:,ii,jj) = sq(phi(:,1,2));
                    out(sc).plph(:,ii,jj) = sq(angle(phloc(:,1,2)));
                    out(sc).ploc(:,ii,jj) = sq(abs(phloc(:,1,2)));
                    out(sc).yerr(:,ii,jj) = sq(yerr(:,1,2,1));
                end
            end
            %close(h);
            out(sc).ElClu = elclu(gfr,:);
            out(sc).FirRate = fr(gfr);
            out(sc).State = States{s};
            sc= sc+1;
        end %loop through states
        save([FileBase '.' mfilename '.mat'],'out');

    case 'display'
        load([FileBase '.' mfilename '.mat']);
        %  global OutArgs;
        %  Map = MapSilicon(x,Channels,ChInShank)

        CluLoc = ClusterLocation(FileBase,El);
        for s=1:length(out)
            fprintf('%s\n',out(s).State);
            if ifTetr==0
                figure(3232);clf
                for jj=1:size(out(1).coh,3)
                    for s=1%:nSt
                        mycell = find(CluLoc(:,1)==out(s).ElClu(jj,1) & CluLoc(:,2)==out(s).ElClu(jj,2));
                        mychan = CluLoc(mycell,3);
                        mychani = find(RepCh==mychan);
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
                        imagesc(out(s).f,[1:nRepCh],sq(out(s).coh(:,:,jj))');hold on
                        colorbar('NorthOutside');
                        title('coherences'); ylabel('channels');
                        hl = Lines([],mychani,'w',[],2);
                        hold off

                        subplot(5,2,5)
                        plot(out(s).f,sq(out(s).powu(:,jj)));axis tight
                        title(['unit power - ' num2str(out(s).ElClu(jj,:)) ' ' ....
                            num2str(out(s).FirRate(jj)) ' Hz, State ' States{s}]);
                        while 1
                            subplot(5,2,[1 3])
                            [x y b] = ginput(1);
                            if b==1
                                [dummy fi] = min(abs(out(s).f-x));
                                %[dummy chi] = min(abs([1:size(out(s).coh,3)]-y));
                                chi = round(y); chi = max(chi,1); chi = min(chi,size(out(s).coh,2));
                                mycoh = sq(out(s).coh(fi,:,jj));
                                subplot(3,2,2)
                                MapSilicon(mycoh,Par);
                                hold on
                                hpt = plot(xch, ych, 'ko','MarkerSize',20);
                                hold off
                                if isfield(out(s), 'ploc')
                                    myploc = sq(out(s).ploc(fi,:,jj));
                                else
                                    myploc = 1;
                                end
                                subplot(3,2,4)
                                MapSilicon(myploc,Par);
                                hold on
                                hpt = plot(xch, ych, 'ko','MarkerSize',20);
                                hold off

                                mypow = sq(out(s).powf(fi,:));
                                %myph = sq(out(s).ph(fi,:,jj));
                                subplot(3,2,6)
                                MapSilicon(mypow,Par);
                                %                            MapSilicon(myph,Par,16);
                                hold on
                                hpt = plot(xch, ych, 'ko','MarkerSize',20);
                                hold off

                                subplot(5,2,7)
                                plot(out(s).f,sq(out(s).coh(:,chi,jj)));axis tight
                                hold on
                                plot(out(s).f,sq(out(s).yerr(:,chi,jj)),'r--');
                                hold off
                                title(['coherence to channel ' num2str(RepCh(chi))]);

                                subplot(5,2,9)
                                plot(out(s).f,sq(out(s).ph(:,chi,jj)));axis tight
                                title(['phase shift to channel ' num2str(RepCh(chi))]);

                            elseif b==3
                                break;
                            else
                                return;
                            end
                        end %while

                    end
                end

            else
                figure(3232);clf
                for jj=1:size(out(s).coh,3)
                    subplotfit(jj,size(out(s).coh,3));
                    imagesc(out(s).f,[1:nRepCh],sq(out(s).coh(:,:,jj))');
                end

            end
            waitforbuttonpress
        end %loop over states
end





