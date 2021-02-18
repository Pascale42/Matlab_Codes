function out = PhaseTimeAnal(FileBase,varargin)
%function out = PhaseTimeAnal(FileBase,State,ElLoc,Display,Overwrite, Report,UseDent)
%does all the phase estimation for all cells from certian location
% UseDent - flag to use dentate phase if it exists
[States, ElLoc,Display, Overwrite, Report,UseDent] = DefaultArgs(varargin,{{'REM','RUN'},'c',1,0,0,0});
out = struct([]);

MinPeriod = 1; %seconds for theta periods selection
MinNumSpks = 100;
WinSec =2;
par = LoadPar([FileBase '.xml']);
sFs = 1e6/par.SampleTime;
eFs = GetEegFs(FileBase);
pFs = 125; % sampling rate of phase vectors - since I did resampling in MakePhaseFiles
warning off
if isstr(ElLoc(1))
    if strcmp(ElLoc,'all')
        El = [1:par.nElecGps];
    else
        El = find(strcmp(par.ElecLoc,ElLoc));
    end
else
    El = ElLoc;
end

if isempty(El)
    fprintf('No electrodes from %s\n',ElLoc);
    return;
end

if ischar(States); States = {States}; end;
nSt = length(States);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%computation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
matname = [FileBase '.pta-' ElLoc ];

if ~FileExists(matname) | Display==2 |Overwrite

    %channels used for phase calculation
    Channels = load([FileBase '.eegseg.par']);
    nCh = length(Channels);
    
        
    %can use hpc and dentate
    if length(Channels)==3
        WithDent = 1; nLoc = 2;
    else
        WithDent = 0; nLoc =1;
    end;

    if ~UseDent;    WithDent=0; end;

    if UseDent & WithDent;   ThCh=3; else     ThCh=1; end

    eeg = readmulti([FileBase '.eeg'],par.nChannels,Channels(ThCh));
    eeg = resample(eeg,1,10);                   %RESAMPLED EEG
    
    %now load all spikes
    [Res, Clu, nClu, spikeph, ClustByEl, cID] = ReadEl4CCG(FileBase,El);
    %keyboard
    clear spikeph; %don't use the old one
    Res =  max(1,round(Res*pFs/sFs)); %get spike times to the sample rate of phase.avoid zero


    %now load the phases vectors (and amplitude)
    ThPh = bload([FileBase '.thphases'], [nCh inf])*pi/32767;
    ThPh = double(ThPh)'; 
    ThTr = LocalMinima(ThPh(:,1),10,0);
    ThPk = LocalMinima(abs(ThPh(:,1)),10,pi/2);
    
    
    %keyboard
    % compute the correction for the phase - later
    % [PhCdf x]= cdfcall(ThPh); bin x; PhaseCorr = 2*pi*PhCdf - pi;
    %get the continuous phase - from hilbert transofrm
    %ContPh = ThPh(Res);
    %ContAmp = ThAmp(Res);

    %compute the phase from troughs - LATER
    %TrPh =

    Out = struct([]);
    for w = 1:length(States) %loop through states

        %load theta periods (more then MinPeriod)
        if FileExists([FileBase '.sts.' States{w}])
            AllPeriod = load([FileBase '.sts.' States{w}]);
            Period = AllPeriod(find(diff(AllPeriod,1,2)>eFs*MinPeriod),:);
            Period = Period*pFs/eFs;
        else
            fprintf('No %s periods. Return empty output\n',States{w});
            continue;
        end
        if isempty(Period);
            fprintf('No %s periods. Return empty output\n',States{w});
            continue;
        end

        PeriodTime = sum(diff(Period,[],2))/pFs; % in seconds
        fprintf('\nFile %s has %2.1f seconds of %s \n',FileBase,PeriodTime,States{w});
 
        %here everything is in spikes sample rate - too slow
%         [out(w).ccg,out(w).t] = Trains2CCG({ThTr*sFs/pFs,Res},{1,Clu},10,30,sFs,'count',Period*sFs/eFs,1);
%         [tr,gr] = CatTrains({ThTr*sFs/pFs,Res},{1,Clu},Period*sFs/eFs,1);
        Period = round(Period);
        Focus = WithinRanges(Res, Period);
        Focus = Focus==1; % make it logical to use fast subscription
        
                      
        
        %select spikes that are in that State
        %         uClu = unique(Clu(Focus));
        %         nClu = length(uClu);
        cnt = 1;
        %         myClu = ismember([1:nClu],uClu); %indexes of clusters that come into cohernce computation (present in this state)
        fprintf('State : %s - PROCESSING\n',States{w});
        
        % fields of struct PhOut :  pR, th0, R, phconfint, mu, k, pNU, V, phhist, phbin
        %now loop through cells
        for c=1:nClu

            %limit Res to the Periods
            myRes = Res(find(Clu==c & Focus));
            if length(myRes)<MinNumSpks
                continue;
            end
            
            myph = ThPh(myRes,ThCh);
            myphc = PhaseFromCycles(myRes,ThPk(1:end-1), ThPk(2:end));
            
            %compute running correlation btw phase and time
            
            [OvSeg, OvSpkInd, OvSegInd, OvSegGrpInd] = FitOverlapedBins(pFs*10,pFs*1,myRes,0);
          %  whos myph
           out(c,w).rc = zeros(size(OvSeg,1),1);
          out(c,w).pv = zeros(size(OvSeg,1),1);
            for ii=1:size(OvSeg,1)
                if sum(OvSegInd==ii)>10
                    spkii = OvSpkInd(OvSegInd==ii);
                    [out(c,w).rc(ii) out(c,w).pv(ii)]= RankCorrelation(myRes(spkii),myph(spkii));
                
                end;
                out(c,w).rct(ii)=OvSeg(ii,1);
            end
            
            
            %compute CCG
            CCGBinSize = round(20*pFs/1000); HalfBins = 30;
            [tr,gr] = CatTrains({ThTr,myRes},{1,2},Period,1);
            %[out(c,w).ccg, out(c,w).t, pairs] = CCG(tr,gr,CCGBinSize,HalfBins,pFs,[1 2],'count');
        
            %[mat, out(c,w).ccgt out(c,w).ccgy]= CCGEvolution(tr, gr, round(40*pFs/1000), 20, pFs, [1 2], WinSec*pFs,Period);
%             [mat, out(c,w).ccgt out(c,w).ccgy]= CCGEvolutionOverlap(tr, gr, CCGBinSize, HalfBins, pFs, [1 2], WinSec*pFs,pFs);
%             out(c,w).ccgev(:,:) = sq(mat(:,:,1,2))';
%             out(c,w).acgev(:,:) = sq(mat(:,:,2,2))';
            out(c,w).cID = cID(c);
            out(c,w).ElClu = ClustByEl(c,:);
            
            % test for any alternative to unofrm
            [H, out(c,w).V, out(c,w).pNU] = TestPhaseNonuniform(myph);
             
            %Step = 1*pFs; Window = 2*pFs;
            %[out(c,w).ph_hist, out(c,w).pht, out(c,w).ph_bins] = RunningHistDiscr(myRes, myph*180/pi, Window,Step , 19);
            %do non-verlapping running hist
            [Starts TimeInd] = FitEvenlySpacedBins(WinSec*pFs, Period, myRes);
            
            nWins = length(Starts); %how many time bins we will have total!
            gi = TimeInd>0; %those that fit in windows
            TimeInd = TimeInd(gi); myph = myph(gi); myphc = myphc(gi);
            
            
            [dummy PhInd] = histcI(myph*180/pi,[-180:20:180]);
            %[d1 d2 NewTimeInd] = unique(TimeInd(gi));
            out(c,w).phhist = Accumulate([TimeInd PhInd],1,[nWins,18]);
            out(c,w).phy = linspace(-180,180,18);
            out(c,w).pht = (Starts'+WinSec/2*pFs)/pFs;
            
            %non-overlapping windows of 2secs
            % this is my way
            %[y, f, t, phi] = mtptchg(eeg, myRes, ones(length(myRes),1), 2^9, pFs, WinSec*pFs, 0 , 3, 'linear', 5, [1 20], [1],Period);
%             out(c,w).coh = sq(y(:,:,1,2));
%             out(c,w).powu = sq(y(:,:,2,2));
%             out(c,w).powf = sq(y(:,:,1,1));
%             out(c,w).phi = sq(phi(:,:,1,2));
%             out(c,w).f = f; out(c,w).t = t;
            %this is using parthas
            [DataCont Starts] = Data2Chronux(eeg, 'c', WinSec*pFs, Period);
            [DataSpk Starts cm]  = Data2Chronux(myRes, 'd', WinSec*pFs, Period, pFs);
            cm = find(cm);
            %get skipped no-spikes windows
            [dummy NewWinInd NewSpkInd] = unique(DataSpk(:,2));
            DataSpk(:,2) = NewSpkInd;
            DataCont = DataCont(:,cm);
            
            %null the matrixes
            
            [C,phi,f,confC,phierr,Cerr]= coherencycpt(DataCont,DataSpk,[3 5],1,pFs,[1 20],[2 0.95] ,0);
            
            sz= [nWins,length(f)];    
            out(c,w).coh = zeros(sz);out(c,w).phi = zeros(sz);out(c,w).cerr = zeros(nWins, length(f),2);out(c,w).phierr=zeros(sz);
            out(c,w).confC = zeros(nWins); out(c,w).powc=zeros(sz);out(c,w).pows=zeros(sz);
            out(c,w).powcerr = zeros(nWins, length(f),2); out(c,w).powserr = zeros(nWins, length(f),2); out(c,w).Fr =zeros(nWins,1);
            
            out(c,w).coh(cm,:) = C';
            out(c,w).phi(cm,:)  = phi'; ,out(c,w).f =f';
            out(c,w).confC(cm) = confC'; 
            out(c,w).cerr(cm,:,:) = permute(Cerr,[3 2 1]); 
            out(c,w).phierr(cm,:)= phierr';
            out(c,w).t = (Starts+WinSec/2*pFs)/pFs;
          
            
            %power of spikes
            [S,f,R,Serr]=mtspectrumpt(DataSpk,[3 5],1,pFs,[1 20],[2 0.95] ,0,1,WinSec*pFs);
            out(c,w).pows(cm,:) = S';
            out(c,w).powserr(cm,:,:) = permute(Serr,[3 2 1]);
            out(c,w).Fr(cm) = R';
            
            
            %power of continuous data (whitened)
            [DataCont Starts] = Data2Chronux(WhitenSignal(eeg), 'c', WinSec*pFs, Period);
            DataCont = DataCont(:,cm);
            [S,f,Serr]=mtspectrumc(DataCont,[3 5],1,pFs,[1 20],[2 0.95],0);
            out(c,w).powc(cm,:) = S';
            out(c,w).powcerr(cm,:,:) = permute(Serr,[3 2 1]);
                                  
                        
            %now compute significance by time
            % test for unimodal vs uniform distr
            %out(c,w).pR=zeros(nWins,1); out(c,w).th0=zeros(nWins,1); out(c,w).R=zeros(nWins,1); out(c,w).logZ=zeros(nWins,1);
            [out(c,w).pR out(c,w).th0 out(c,w).R out(c,w).logZ] = RayleighTest(myph,TimeInd);
            out(c,w).th0 = out(c,w).th0*180/pi;
              curl = length(out(c,w).pR);
                 if curl <nWins
                    out(c,w).pR(curl+1:nWins)=zeros(nWins-curl,1); 
                    out(c,w).th0(curl+1:nWins)=zeros(nWins-curl,1); 
                    out(c,w).logZ(curl+1:nWins)=zeros(nWins-curl,1); 
                    out(c,w).R(curl+1:nWins)=zeros(nWins-curl,1); 
                 end
            
            [out(c,w).pRc out(c,w).th0c out(c,w).Rc out(c,w).logZc] = RayleighTest(myphc,TimeInd);
             out(c,w).th0c = out(c,w).th0c*180/pi;
              curl = length(out(c,w).pRc);
                 if curl <nWins
                    out(c,w).pRc(curl+1:nWins)=zeros(nWins-curl,1); 
                    out(c,w).th0c(curl+1:nWins)=zeros(nWins-curl,1); 
                    out(c,w).logZc(curl+1:nWins)=zeros(nWins-curl,1); 
                    out(c,w).Rc(curl+1:nWins)=zeros(nWins-curl,1); 
                 end
             
            [PhSegs StartsSeg] = Data2Chronux(ThPh(:,ThCh), 'c', WinSec*pFs, Period);
            PhSegs= PhSegs(:,cm);
            sc =linspace(0.8, 1.2, 19);%[3:0.3:15]/8;
            [SpkSegs SpkStarts]  = Data2Chronux(myRes, 'd', WinSec*pFs, Period, 1);
            [dummy NewWinInd NewSpkInd] = unique(SpkSegs(:,2));
            SpkSegs(:,2) = NewSpkInd;
            
            %[Starts TimeInd] = FitEvenlySpacedBins(WinSec*pFs, Period, myRes);
            %PhSegs = GetSegs(ThPh(:,ThCh),Starts,WinSeg*pFs
            for ii=1:size(PhSegs,2) % loop through segments
                if sum(SpkSegs(:,2)==ii)>1
                    PhResc = PhaseRescale(PhSegs(:,ii), SpkSegs(SpkSegs(:,2)==ii,1), sc,0,0);
                    [out(c,w).pRsc(cm(ii),:) out(c,w).th0sc(cm(ii),:) out(c,w).Rsc(cm(ii),:) out(c,w).logZsc(cm(ii),:)] = RayleighTest(PhResc);              
                else
                    out(c,w).pRsc(cm(ii),:) =ones(length(sc),1)';
                    out(c,w).th0sc(cm(ii),:) =zeros(length(sc),1)';
                    out(c,w).Rsc(cm(ii),:) =zeros(length(sc),1)';
                    out(c,w).logZsc(cm(ii),:) = zeros(length(sc),1)';
                end
                    
            end
            out(c,w).scy = sc;
            out(c,w).sct = StartsSeg(cm);
            
                       
                fprintf('State: %s Cell %10f   \n',...
                States{w}, out(c,w).cID );

            
        end

    end
    if Overwrite
        save(matname,'out');
    end
else
    
    load(matname,'-MAT');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%display
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if Display>0
    [nClu, nSt] = size(out);
    
    if 0
    
    for c=1:nClu
        anyc=0;
        for w = 1:nSt %loop through states
            if isempty(out(c,w).cID);  continue;    end
    %        anyc =anyc+1;
            
            %compute the freq. of each segment accoring to rescaling
            gsc=find(out(c,w).scy>0.7 & out(c,w).scy<1.3);
            [maxr msc] = max(out(c,w).Rsc(:,gsc),[],2);
            ppst = PowerPeakStats(out(c,w).powc(1:length(msc),:),out(c,w).f,[5 12]);
            maxsc = out(c,w).scy(gsc(msc));
            thfr = sq(ppst(1,:,1));
            scfr = thfr.*maxsc;
                        
            figure(876)
            clf
        
            fprintf('State: %s Cell: El=%d, Clu=%d,\n', States{w}, out(c,w).ElClu(1), out(c,w).ElClu(2));
            npl=5;
            subplot(3, 2, 1)
            smimage(out(c,w).ccgt,out(c,w).ccgy,sq(out(c,w).ccgev),2,0);
            
            subplot(3, 2, 2)
            smimage(out(c,w).ccgt,out(c,w).ccgy,sq(out(c,w).acgev),2,0);
                        
            subplot(3, 2, 3)
            %smimage(out(c,w).pht,out(c,w).phy,out(c,w).phhist,2,0.7);
            smimage(out(c,w).t,out(c,w).scy,unity(out(c,w).Rsc')',2,0.9);
            
            subplot(3, 2, 4)
            scatter(out(c,w).t(1:length(msc)),scfr,5,maxr);axis tight
            %xlim([min(out(c,w).t) max(out(c,w).t)]);
            
            subplot(3, 2, 5)
            smimage(out(c,w).t,out(c,w).f,abs(out(c,w).pows),2,0.8);axis xy
            
            subplot(3, 2, 6)
            smimage(out(c,w).t,out(c,w).f,out(c,w).powc,2,0.8);axis xy;ylim([4 13]);
            %smimage(out(c,w).t,out(c,w).f,out(c,w).coh,2,0.8);
            %plot phase hist
            %subplot(nSt,4,1+(w-1)*4)
            
            drawnow
            TimeBrowse(100,100); 
        
                               
            %keyboard
        
            %waitforbuttonpress
                   

        end
        
        
%         if anyc>0%out(c,w).logZ>1.1 || out(c,w).V>1.7 % <0.0001 %| out(c,w).k>0.3
% 
%             while 1
%                 waitforbuttonpress
%                 chld = get(gcf,'Children');
%                 %if nSt==1
%                 chld = flipud(chld(:));
%                 %else
%                 %   chld = reshape(flipud(reshape(chld(:),[],2)),[],1);
%                 %end
%                 sel = get(gcf,'SelectionType');
%                 mousecoord = get(gca,'CurrentPoint');
%                 xmouse = mousecoord(1,1);
%                 ymouse = mousecoord(1,2);
%                 ca = get(gcf,'CurrentAxes');
%                 curwin = find(chld==ca);
%                 if curwin>6; curwin = curwin-6; end
% 
%                 switch sel
%                     case 'normal'
%                         break
%                     case 'alt'
%                         for curst=1:nSt
%                             if isempty(out(c,curst).cID);  continue;    end
%                             subplot(nSt,4,1+(curst-1)*4)
%                             if curwin==2 || curwin==3
%                                 [dummy curlag] = min(abs(out(c,curst).sh_lag-xmouse));
%                                 bar([out(c,curst).phbin out(c,curst).phbin+360],[out(c,curst).sh_hist(:,curlag); out(c,curst).sh_hist(:,curlag)]);
%                                 axis tight
%                             elseif curwin==5
%                                 [dummy curfr] = min(abs(out(c,curst).fr_yb(1:end-1)-ymouse));
%                                 bar(out(c,curst).fr_xb(1:end-1), out(c,curst).fr_hist(:,curfr));
%                                 axis tight
% 
%                             elseif curwin==6
%                                 [dummy curamp] = min(abs(out(c,curst).amp_yb(1:end-1)-ymouse));
%                                 bar(out(c,curst).amp_xb(1:end-1), out(c,curst).amp_hist(:,curamp));
%                                 axis tight
%                             end
%                         end
%                         %keyboard
%                 end
%             end
% 
%             if Report
%                 rstr = sprintf('Cell %10f : log(Z)=%2.2f pR=%2.4f, R=%2.2f, k=%2.2f, pNU=%2.2f , V=%2.2f\n',...
%                     out(c,w).cID, out(c,w).logZ, out(c,w).pR, out(c,w).R, out(c,w).k, out(c,w).pNU, out(c,w).V);
% 
%                 reportfig(gcf, 'PhaseAnal',0,[FileBase '  ' rstr ],150);
%             end
%         end

    end
    end %if 0
    anyc=0;
    for w = 1:nSt %loop through states
        figure(212+w);clf; figure(213+w);clf
        set(212+w,'RendererMode','manual');
        set(213+w,'RendererMode','manual');
        set(212+w,'Renderer','zbuffer');
        set(213+w,'Renderer','zbuffer');
        cm = colormap;
        ci = 1+[0:nClu-1]*round(size(cm,1)/nClu);
        for c=1:nClu
            
          if isempty(out(c,w).cID) | isempty(out(c,w).ccgev);  continue;    end
          %figure(212+w)
          % subplotfit(c,nClu);
           %[result(c,w).maxlag result(c,w).tt] = 
%           subplot(311)
           
           
           PeaksRidge(out(c,w).ccgev,[2 2],20,out(c,w).ccgt,out(c,w).ccgy, [0 0 0 0 0]); 
           tstr = sprintf('%s: El=%d,Clu=%d', States{w}, out(c,w).ElClu(1), out(c,w).ElClu(2));
           %pause
           title(tstr);
           waitforbuttonpress;
           if 0
               figure(213+w)
               subplot(211)
               plot(result(c,w).tt,result(c,w).maxlag,'-*','Color',cm(ci(c),:));hold on
               subplot(212)
               plot(result(c,w).tt(1:end-1),abs(diff(result(c,w).maxlag)),'-*','Color',cm(ci(c),:));hold on
           end
%            subplot(223)
%            plot(result(c,w).tt,result(c,w).maxlag,'-*','Color',cm(ci(c),:));
        end
        
    %$pause
    
    %figure(212)
%     reportfig(212,'RidgeAnal',0,FileBase,200);
%     reportfig(213,'RidgeAnal',0,FileBase,200);
     end
end

%Out = CatStruct(out,[],3);





