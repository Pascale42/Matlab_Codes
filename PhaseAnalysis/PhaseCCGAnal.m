function out = PhaseCCGAnal(FileBase,varargin)
%function out = PhaseCCGAnal(FileBase,State,ElLoc,Display,Overwrite, Report,UseDent)
%does all the phase estimation for all cells from certian location
% UseDent - flag to use dentate phase if it exists
[States, ElLoc,Display, Overwrite, Report,UseDent] = DefaultArgs(varargin,{'REM','c',1,0,0,0});
out = struct([]);

MinPeriod = 1; %seconds for theta periods selection
MinNumSpks = 100;

par = LoadPar([FileBase '.par']);
sFs = 1e6/par.SampleTime;
eFs = GetEegFs(FileBase);
pFs = 125; % sampling rate of phase vectors - since I did resampling in MakePhaseFiles
warning off
if isstr(ElLoc(1))
    El = find(strcmp(par.ElecLoc,ElLoc));
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
if ~FileExists([FileBase '.pccg']) | Overwrite

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

    %now load all spikes
    [Res, Clu, nClu, spikeph, ClustByEl, cID] = ReadEl4CCG(FileBase,El);
    %keyboard
    clear spikeph; %don't use the old one
    Res =  max(1,round(Res*pFs/sFs)); %get spike times to the sample rate of phase.avoid zero


    %now load the phases vectors (and amplitude)
    ThPh = bload([FileBase '.thphases'], [nCh inf])*pi/32767;
    ThAmp = bload([FileBase '.thamps'], [nCh inf])/32767;
    ThPh = ThPh'; ThAmp = ThAmp';
    k = 4; % 2
    gaussker = normpdf(-k:k,0, k/4)';
    sThPh = conv2(gaussker,1,unwrap(ThPh),'same');
    ThFr = diff(sThPh)*pFs/2/pi;

    ThTr = LocalMinima(ThPh(:,1),10,0);

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
            
        
%             [mat, out(w).evxr out(w).evyr]= CCGEvolution(tr, gr, 20*sFs/1000, 20, sFs, [1 c+1], 2*sFs);
%             out(w).ccgev(:,:,c) = sq(mat(:,:,1,2));

            out(c,w).cID = cID(c);
            out(c,w).ElClu = ClustByEl(c,:);
            myph = ThPh(myRes,ThCh);
            myfr = ThFr(myRes,ThCh);
            myamp = log(ThAmp(myRes,ThCh));
            % test for unimodal vs uniform distr
            [out(c,w).pR out(c,w).th0 out(c,w).R out(c,w).logZ] = RayleighTest(myph);

            % test for any alternative to unofrm
            [H, out(c,w).V, out(c,w).pNU] = TestPhaseNonuniform(myph);

            %out(c,w).logZ = log(-log(out(c,w).pval));
            out(c,w).th0 = out(c,w).th0*180/pi;

            [cm, r] = CircConf(myph,0.01,200);
            out(c,w).phconfint = r*180/pi;
            [out(c,w).mu out(c,w).k] = VonMisesFit(myph);

            out(c,w).n = length(myRes);

            %cvPh = PhaseFit(myph);
            %out(c,w).cvL1 = cvPh.L1;
            %out(c,w).cvL2 = cvPh.L2;
            %myph = myph*180/pi;
            %myphall = mod(myphall*180/pi,360);

            edges = [-180:20:180];
            [out(c,w).phhist ] = histcI(myph*180/pi,edges);
            out(c,w).phbin = (edges(2:end)+edges(1:end-1))/2;
            out(c,w).phhist = out(c,w).phhist';

            %compute lagged phase modulation

            MaxShift = round(pFs*2000/1000); BinSize = round(pFs/1000*20);
            [out(c,w).sh_p, out(c,w).sh_lag, out(c,w).sh_th0, out(c,w).sh_r, out(c,w).sh_logZ, out(c,w).sh_hist] = ...
                ShiftPhase(ThPh(:,ThCh),myRes,MaxShift,BinSize,Period);
            out(c,w).sh_lag =out(c,w).sh_lag/pFs*1000;

            %compute shifted freq.
            [out(c,w).sh_frav,out(c,w).sh_frlag, out(c,w).sh_frstd] = ShiftParam(ThFr(:,ThCh),myRes,MaxShift,BinSize,Period);
            % .... and amplitude of theta
            [out(c,w).sh_ampav,out(c,w).sh_amplag, out(c,w).sh_ampstd] = ShiftParam(ThAmp(:,ThCh),myRes,MaxShift,BinSize,Period);

                                   
            %compute phase histogram as a function of freq.
            prct1 = prctile(myfr,[5 95]); gi = find(myfr>prct1(1) & myfr<prct1(2));
            [out(c,w).fr_hist, out(c,w).fr_xb,out(c,w).fr_yb] = hist2([myph(gi)*180/pi myfr(gi); myph(gi)*180/pi+360 myfr(gi)], 18, 15);
            
            %compute circ regression 
            [out(c,w).crd_fr out(c,w).cru_fr]= CircCorr(myph(gi), myfr(gi));
                        
            %compute phase histogram as a function of amp.
            prct2 = prctile(myamp,[5 95]); gi = find(myamp>prct2(1) & myamp<prct2(2));
            [out(c,w).amp_hist, out(c,w).amp_xb,out(c,w).amp_yb] = hist2([myph(gi)*180/pi myamp(gi); myph(gi)*180/pi+360 myamp(gi)], 18, 10);

            %compute circ regression 
            [out(c,w).crd_amp out(c,w).cru_amp]= CircCorr(myph(gi), myamp(gi));
       
            %compute CCG
            CCGBinSize = round(20*pFs/1000); HalfBins = 30;
            [tr,gr] = CatTrains({ThTr,myRes},{1,2},Period,1);
            [out(c,w).ccg, out(c,w).t, pairs] = CCG(tr,gr,CCGBinSize,HalfBins,pFs,[1 2],'count');
            
            
            %now compute the CCG sorted by frequency
            MyPairs = pairs(gr(pairs(:,1))==1 & gr(pairs(:,2))==2,:);
            BinInd = round(diff(tr(MyPairs),1,2)/CCGBinSize) + HalfBins + 1;
            
            nbins = 10; %for amp and freq.
            bins = linspace(prct1(1),prct1(2),nbins+1);
            [FrHist FrInd] = histcI(ThFr(tr(MyPairs(:,2))), bins);
            gbins = find(FrInd>0);
            CountFr = Accumulate([BinInd(gbins), FrInd(gbins)], 1, [2*HalfBins+1, nbins]);
            
            Eps =1;
            % scale it with Poisson assumption
            n1 = histcI(ThFr(tr(gr==1)), bins);
            n2 = histcI(ThFr(tr(gr==2)), bins);
            Expected = n1.*n2/(2*HalfBins+1);
            out(c,w).ccg_fr = CountFr;%(CountFr+Eps) ./(repmat(Expected',2*HalfBins+1,1)+Eps);
            
            %same for the amplitude
            bins = linspace(prct2(1),prct2(2),nbins+1);
            [AmpHist AmpInd] = histcI(log(ThAmp(tr(MyPairs(:,2)),ThCh)), bins);
            gbins = find(AmpInd>0);
            CountAmp= Accumulate([BinInd(gbins), AmpInd(gbins)], 1, [2*HalfBins+1, nbins]);
            
            n1 = histcI(ThAmp(tr(gr==1)), bins);
            n2 = histcI(ThAmp(tr(gr==2)), bins);
            Expected = n1.*n2/(2*HalfBins+1);
            out(c,w).ccg_amp = CountAmp; %(CountAmp+Eps) ./(repmat(Expected',2*HalfBins+1,1)+Eps);
            
            
            fprintf('State: %s Cell %10f : n=%d, log(Z)=%2.2f pR=%2.4f, R=%2.2f, k=%2.2f, pNU=%2.2f , V=%2.2f\n',...
                States{w}, out(c,w).cID, out(c,w).n, out(c,w).logZ, ...
                out(c,w).pR, out(c,w).R, out(c,w).k, out(c,w).pNU, out(c,w).V);


            %Out(w) = a;
            %a=CatStruct(out);
        end

    end
    save([FileBase '.pccg'],'out');
else
    load([FileBase '.pccg'],'-MAT');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%display
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if Display
    [nClu, nSt] = size(out);

    for c=1:nClu
        figure(876)
        clf
        anyc=0;
        for w = 1:nSt %loop through states
            if isempty(out(c,w).cID);  continue;    end
            anyc =anyc+1;
            fprintf('State: %s Cell: El=%d, Clu=%d, n=%d, log(Z)=%2.2f pR=%2.4f, R=%2.2f, k=%2.2f, pNU=%2.2f , V=%2.2f\n',...
                States{w}, out(c,w).ElClu(1), out(c,w).ElClu(2), out(c,w).n, out(c,w).logZ, out(c,w).pR, out(c,w).R, out(c,w).k, out(c,w).pNU, out(c,w).V);


            %plot phase hist
            subplot(nSt,4,1+(w-1)*4)
            bar([out(c,w).phbin out(c,w).phbin+360],[out(c,w).phhist out(c,w).phhist]);
            axis tight
            ylabel(States{w});

            %plot shift modulation
            %subplot(nSt,4,2+(w-1)*4)
            subplot(nSt*2,4,2+(w-1)*8)
            plot(out(c,w).sh_lag,out(c,w).sh_logZ);
            axis tight
            yl = ylim; ylim([-1 yl(2)]);

            %plot shift params
            subplot(nSt*2,4,6+(w-1)*8)
            ax = plotyy(out(c,w).sh_frlag/pFs*1000,out(c,w).sh_frav,out(c,w).sh_amplag/pFs*1000,out(c,w).sh_ampav);
            hold(ax(1),'on'); hold(ax(2),'on');
            %plot(ax(1), out(c,w).sh_frlag,out(c,w).sh_frav+out(c,w).sh_frstd,'--');
            %plot(ax(2), out(c,w).sh_amplag,out(c,w).sh_ampav+out(c,w).sh_ampstd,'--');
            axis(ax(1),'tight');
            axis(ax(2),'tight');

            %plot fr. hist
            subplot(nSt,4,3+(w-1)*4)
            imagesc(out(c,w).fr_xb(1:end-1), out(c,w).fr_yb(1:end-1),unity(out(c,w).fr_hist)');
            axis tight; axis xy

            %plot time
            subplot(nSt,4,4+(w-1)*4)
            imagesc(out(c,w).amp_xb(1:end-1), out(c,w).amp_yb(1:end-1),unity(out(c,w).amp_hist)');
            axis tight; axis xy

        end
        drawnow
        if anyc>0%out(c,w).logZ>1.1 || out(c,w).V>1.7 % <0.0001 %| out(c,w).k>0.3

            while 1
                waitforbuttonpress
                chld = get(gcf,'Children');
                %if nSt==1
                chld = flipud(chld(:));
                %else
                %   chld = reshape(flipud(reshape(chld(:),[],2)),[],1);
                %end
                sel = get(gcf,'SelectionType');
                mousecoord = get(gca,'CurrentPoint');
                xmouse = mousecoord(1,1);
                ymouse = mousecoord(1,2);
                ca = get(gcf,'CurrentAxes');
                curwin = find(chld==ca);
                if curwin>6; curwin = curwin-6; end

                switch sel
                    case 'normal'
                        break
                    case 'alt'
                        for curst=1:nSt
                            if isempty(out(c,curst).cID);  continue;    end
                            subplot(nSt,4,1+(curst-1)*4)
                            if curwin==2 || curwin==3
                                [dummy curlag] = min(abs(out(c,curst).sh_lag-xmouse));
                                bar([out(c,curst).phbin out(c,curst).phbin+360],[out(c,curst).sh_hist(:,curlag); out(c,curst).sh_hist(:,curlag)]);
                                axis tight
                            elseif curwin==5
                                [dummy curfr] = min(abs(out(c,curst).fr_yb(1:end-1)-ymouse));
                                bar(out(c,curst).fr_xb(1:end-1), out(c,curst).fr_hist(:,curfr));
                                axis tight

                            elseif curwin==6
                                [dummy curamp] = min(abs(out(c,curst).amp_yb(1:end-1)-ymouse));
                                bar(out(c,curst).amp_xb(1:end-1), out(c,curst).amp_hist(:,curamp));
                                axis tight
                            end
                        end
                        %keyboard
                end
            end

            if Report
                rstr = sprintf('Cell %10f : log(Z)=%2.2f pR=%2.4f, R=%2.2f, k=%2.2f, pNU=%2.2f , V=%2.2f\n',...
                    out(c,w).cID, out(c,w).logZ, out(c,w).pR, out(c,w).R, out(c,w).k, out(c,w).pNU, out(c,w).V);

                reportfig(gcf, 'PhaseAnal',0,[FileBase '  ' rstr ],150);
            end
        end

    end
end

%Out = CatStruct(out,[],3);





