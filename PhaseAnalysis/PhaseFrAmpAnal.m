function out = PhaseFrAmpAnal(FileBase,varargin)
%function out = PhaseFrAmpAnal(FileBase,fMode,State,ElLoc,Overwrite, Report,UseDent)
%does all the phase estimation for all cells from certian location
%Latest version March12 2007
[fMode, States, ElLoc, Overwrite, Report, UseDent] = DefaultArgs(varargin,{'display',{'REM','RUN'},'c',0,0,0});

smx = 1; % length of smoother in x dimmension of histograms for display

Par = LoadPar([FileBase '.xml']);
MinPeriod = 1; %seconds for theta periods selection
MinNumSpks = 50;
%ThFreqRange = [5 12];
sFs = 1e6/Par.SampleTime;
eFs = GetEegFs(FileBase);
%pFs = 125; % sampling rate of phase vectors - since I did resampling in MakePhaseFiles
warning off
if isstr(ElLoc(1))
    El = find(strcmp(Par.ElecLoc,ElLoc));
else
    El = ElLoc;
end

if isempty(El)
    fprintf('No electrodes from %s\n',ElLoc);
    return;
end

if ischar(States); States = {States}; end;
nSt = length(States);

switch fMode
    case {'compute','collect'}
        out = struct([]);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %computation
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %  if ~FileExists([FileBase '.pframp-' ElLoc])

        %channels used for phase calculation
       % Channels = load([FileBase '.eegseg.par'])+1;
  
       %nCh = length(Channels);

        %now load all spikes
        [Res, Clu, Map] = LoadCluRes(FileBase,El);
        %            [Res, Clu, nClu, spikeph, ClustByEl, cID] = ReadEl4CCG(FileBase,El);
        nClu = max(Clu);
        %keyboard
        clear spikeph; %don't use the old one
        Res =  max(1,round(Res*eFs/sFs)); %get spike times to the sample rate of phase.avoid zero

        [ThPhAll, ThAmpAll, ThFrAll] = ThetaParams(FileBase);


        out = struct([]);
        cnt=1;
        for w = 1:length(States) %loop through states

            %load theta periods (all that are longer than MinPeriod)
            if FileExists([FileBase '.sts.' States{w}])
                AllPeriod = load([FileBase '.sts.' States{w}]);
                Period = AllPeriod(find(diff(AllPeriod,1,2)>eFs*MinPeriod),:);
                Period = Period*eFs/eFs;
            else
                fprintf('No %s periods. Return empty output\n',States{w});
                continue;
            end
            if isempty(Period);
                fprintf('No %s periods. Return empty output\n',States{w});
                continue;
            end

            PeriodTime = sum(diff(Period,[],2))/eFs; % in seconds
            fprintf('\nFile %s has %2.1f seconds of %s \n',FileBase,PeriodTime,States{w});

            ThPh = SelectPeriods(ThPhAll,Period,'c',1);
            ThFr = SelectPeriods(ThFrAll,Period,'c',1);
            ThAmp = SelectPeriods(ThAmpAll,Period,'c',1);

            ThPhu = MakeUniformDistr(ThPh); % makes phase uniform in case it isn't

            fprintf('State : %s - PROCESSING\n',States{w});

            % fields of struct PhOut :  pR, th0, R, phconfint, mu, k, pNU, V, phhist, phbin
            %now loop through cells
            for c=1:nClu

                %limit Res to the Periods
                myRes = SelectPeriods(Res(Clu==c), Period,'d',1,1);
                if length(myRes)<=MinNumSpks
                    continue;
                end

                out(cnt).FileBase = FileBase;
                out(cnt).El = Map(c,2);
                out(cnt).Clu = Map(c,3);
                out(cnt).State = States{w};
                out(cnt).StateTime = PeriodTime;

                myph = ThPh(myRes);
                myphu = ThPhu(myRes);
                myphd = myph*180/pi;
                myfr = ThFr(myRes);
                myamp = ThAmp(myRes);

                nspk = length(myph);

                % test for unimodal vs uniform distr
                [out(cnt).pR out(cnt).th0 out(cnt).R out(cnt).logZ ] = RayleighTest(myph);
                [out(cnt).pRu out(cnt).th0u out(cnt).Ru out(cnt).logZu ] = RayleighTest(myphu);
                %            [out(cnt).RT ] = RayleighTest(myph);
                %           [out(cnt).RTu ] = RayleighTest(myphu);
                
                if strcmp(fMode,'collect')
                    out(cnt).ph = {myph(:)};
                    cnt=cnt+1;
                    continue;
                end
                
                
                %now compute the same thing for subsets of spikes
                for k=1:500
                    phsmpl(:,k) = myph(randsample(nspk, 40));
                end
                out(cnt).rtest = RayleighTest(phsmpl);

                % out(cnt).rtest = %DistrStruct(rtest);

                % test for any alternative to unofrm
                [H, out(cnt).V, out(cnt).pNU] = TestPhaseNonuniform(myph);

                %out(cnt).logZ = log(-log(out(cnt).pval));
                out(cnt).th0 = out(cnt).th0*180/pi;

                [cm, r] = CircConf(myph,0.01,200);
                out(cnt).phconfint = r*180/pi;
                [out(cnt).mu out(cnt).k] = VonMisesFit(myph);

                out(cnt).n = length(myRes);

                edges = [-180:20:180];
                [out(cnt).phhist ] = histcI(myph*180/pi,edges);
                out(cnt).phbin = (edges(2:end)+edges(1:end-1))/2;
                out(cnt).phhist = out(cnt).phhist';

                %compute lagged phase modulation - in fact it is
                %testing time lagging of the spikes relative to the
                %field - same as linear relationship btw phase shift and
                %frequency in spectral analysis

                MaxShift = round(eFs*2000/1000); BinSize = round(eFs/1000*20);

                [out(cnt).sh_p, out(cnt).sh_lag, out(cnt).sh_th0, out(cnt).sh_r, out(cnt).sh_logZ, out(cnt).sh_hist] = ...
                    ShiftPhase(ThPh,myRes,MaxShift,BinSize);
                out(cnt).sh_lag =out(cnt).sh_lag/eFs*1000;

                %compute shifted freq.
                [out(cnt).sh_frav,out(cnt).sh_frlag, out(cnt).sh_frstd] = ShiftParam(ThFr,myRes,MaxShift,BinSize);
                % .... and amplitude of theta
                [out(cnt).sh_ampav,out(cnt).sh_amplag, out(cnt).sh_ampstd] = ShiftParam(ThAmp,myRes,MaxShift,BinSize);


                %  compute the freq-amp preference of the cell
                [out(cnt).framp_hist out(cnt).framp_xb out(cnt).framp_yb] = hist2([myfr myamp], 10, 10);


                %compute phase histogram as a function of freq.
                %prct1 = prctile(myfr,[5 95]); gi = find(myfr>prct1(1) & myfr<prct1(2));
                [out(cnt).fr_hist, out(cnt).fr_xb,out(cnt).fr_yb] = hist2([myphd myfr; myphd+360 myfr], 18, 10);

                %compute circ regression
                % [out(cnt).crd_fr out(cnt).cru_fr]= CircCorr(myph, myfr);


                %compute phase histogram as a function of freq.
                %prct2 = prctile(myamp,[5 95]); gi = find(myamp>prct2(1) & myamp<prct2(2));
                [out(cnt).amp_hist, out(cnt).amp_xb,out(cnt).amp_yb] = hist2([myphd myamp; myphd+360 myamp], 18, 10);

                %compute circ regression
                [out(cnt).crd_amp out(cnt).cru_amp]= CircCorr(myph, myamp);

                %                     %compute CCG
                %                     CCGBinSize = round(10*eFs/1000); HalfBins = 60;
                %                     [tr,gr] = CatTrains({ThTr,myRes},{1,2});
                %                     [out(cnt).ccg, out(cnt).t, pairs] = CCG(tr,gr,CCGBinSize,HalfBins,eFs,[1 2],'count');
                %
                %
                %                     %now compute the CCG sorted by frequency
                %                     MyPairs = pairs(gr(pairs(:,1))==1 & gr(pairs(:,2))==2,:);
                %                     BinInd = round(diff(tr(MyPairs),1,2)/CCGBinSize) + HalfBins + 1;
                %
                nbins = 10; %for  freq.
                frbins = linspace(min(ThFr),max(ThFr),nbins+1);
                %                     [FrHist FrInd] = histcI(ThFr(tr(MyPairs(:,2))), frbins);
                %                     gbins = find(FrInd>0);
                %                     CountFr = Accumulate([BinInd(gbins), FrInd(gbins)], 1, [2*HalfBins+1, nbins]);
                %
                %                     Eps =1;
                %                     % scale it with Poisson assumption
                %                     n1 = histcI(ThFr(tr(gr==1)), frbins);
                %                     n2 = histcI(ThFr(tr(gr==2)), frbins);
                %                     Expected = n1.*n2/(2*HalfBins+1);
                %                     out(cnt).ccg_fr = CountFr;%(CountFr+Eps) ./(repmat(Expected',2*HalfBins+1,1)+Eps);
                %
                %                     %same for the amplitude
                %                     %bins = linspace(prct2(1),prct2(2),nbins+1);
                ampbins = linspace(min(ThAmp),max(ThAmp),nbins+1);
                %                     [AmpHist AmpInd] = histcI(ThAmp(tr(MyPairs(:,2))), ampbins);
                %                     gbins = find(AmpInd>0);
                %                     CountAmp= Accumulate([BinInd(gbins), AmpInd(gbins)], 1, [2*HalfBins+1, nbins]);
                %
                %                     n1 = histcI(ThAmp(tr(gr==1)), ampbins);
                %                     n2 = histcI(ThAmp(tr(gr==2)), ampbins);
                %                     Expected = n1.*n2/(2*HalfBins+1);
                %                     out(cnt).ccg_amp = CountAmp; %(CountAmp+Eps) ./(repmat(Expected',2*HalfBins+1,1)+Eps);


                %now compute the resultant vector in the small overlapping windows , as
                %well as freq and amplitude of theta for more detailed analysis
                %later
                if 0
                    BinSize = eFs*6; Step = 3*eFs;
                    [Segs, SpkInd, SegInd] = FitOverlapedBins(BinSize,Step,myRes,1,length(ThPh));

                    %init the dyn struct array
                    out(cnt).dyn = InitStructArray({'ifr','n','t','p','th0','r','logZ','fr','amp','k'},zeros(max(SegInd),1)*NaN);

                    %count spikes per bin
                    out(cnt).dyn.n = Accumulate(SegInd, 1, size(Segs,1));
                    out(cnt).dyn.ifr = out(cnt).dyn.n./BinSize*eFs;

                    goodi = find(out(cnt).dyn.n>1);
                    out(cnt).dyn.t = mean(Segs,2);

                    out(cnt).dyn = CopyStruct(RayleighTest(myph(SpkInd),SegInd),out(cnt).dyn);

                    segfr = Accumulate(SegInd, myfr(SpkInd), size(Segs,1));channels
                    segamp = Accumulate(SegInd, myamp(SpkInd), size(Segs,1));
                    out(cnt).dyn.fr(goodi) = (segfr(goodi)./out(cnt).dyn.n(goodi))';
                    out(cnt).dyn.amp(goodi) = (segamp(goodi)./out(cnt).dyn.n(goodi))';



                    %now compute rayleigh test on the spikes grouped according to
                    %the significance of the modulation

                    %make evenly populated binning of significance (logZ) - 10 bins
                    %            logZ_rg  = prctile(out(cnt).dyn.logZ(goodi),linspace(0,100,11));
                    logZ_rg  =  prctile(log(out(cnt).dyn.r(goodi)),linspace(0,100,11));
                    [seg_num seg_prct] = histcI(log(out(cnt).dyn.r), logZ_rg);
                    out(cnt).phbysig = RayleighTest(myph(SpkInd),seg_prct(SegInd));
                    out(cnt).phbysig.bins =((logZ_rg(1:end-1)+logZ_rg(2:end))/2)';

                    %do the same thing with spikes assigned at random to the "bins"
                    %.. no temporal contiguity
                    tmpph = myph(SpkInd);
                    sigr=struct([]);
                    mfds = {'p','th0','r','logZ','k'};
                    for kk=1:200
                        rsample = randperm(length(tmpph));
                        dynr = RayleighTest(tmpph(rsample),SegInd);
                        [seg_num seg_prct] = histcI(log(dynr.r), logZ_rg);
                        sigr(kk).rt = InitStructArray(mfds,zeros(10,1)*NaN);
                        sigr(kk).rt = CopyStruct(RayleighTest(tmpph(rsample),seg_prct(SegInd)),sigr(kk).rt);

                    end

                    tmp =CatStruct(sigr);
                    sigr_st = CatStruct(tmp.rt,[],2);
                    out(cnt).sigr_m = FunOfStruct(sigr_st,'nanmean',mfds,2);
                    out(cnt).sigr_std = FunOfStruct(FunOfStruct(FunOfStruct(sigr_st,'transpose'),'nanstd',mfds),'transpose');

                    %compute emptrical p-values form random distribution
                    for f=1:length(mfds)
                        randval = sigr_st.(mfds{f});
                        empval = out(cnt).phbysig.(mfds{f});
                        [s si] = sort([empval'; randval'],1);
                        [nbelow dummy] = find(si==1);
                        out(cnt).phbysig_pval.(mfds{f})=1-nbelow/size(s,1);
                    end
                end
                %now test significance of modulation as function of freq. and
                %amp
                out(cnt).phfr = InitStructArray({'p','th0','r','logZ','n','bin','k'}, zeros(length(frbins)-1,1)*NaN);
                out(cnt).phfr.bin = ((frbins(1:end-1)+frbins(2:end))/2)';
                [out(cnt).phfr.n ind] = histcI(myfr,frbins);
                gi = find(ind>0);
                mi = max(ind(gi));
                [out(cnt).phfr.p(1:mi), out(cnt).phfr.th0(1:mi), out(cnt).phfr.r(1:mi), out(cnt).phfr.logZ(1:mi), out(cnt).phfr.k(1:mi) ] ...
                    = RayleighTest(myph(gi),ind(gi));

                out(cnt).phamp = InitStructArray({'p','th0','r','logZ','n','bin','k'},length(ampbins)-1);
                out(cnt).phamp.bin = ((ampbins(1:end-1)+ampbins(2:end))/2)';
                [out(cnt).phamp.n ind] = histcI(myamp,ampbins);
                gi = find(ind>0);
                mi = max(ind(gi));
                [out(cnt).phamp.p(1:mi), out(cnt).phamp.th0(1:mi), out(cnt).phamp.r(1:mi), out(cnt).phamp.logZ(1:mi),out(cnt).phamp.k(1:mi)] ...
                    = RayleighTest(myph(gi),ind(gi));



                fprintf('State: %s ',States{w});
                PrintStruct(out(cnt),[],2);
                cnt=cnt+1;

            end
        end
        if Overwrite
            save([FileBase '.pframp-' ElLoc],'out');%,'ThFrSpec','ThAmpSpec');
        elseif nargout==0

            fprintf('you are wasting cpu, do you want to save?.. saving in tmpfile \n');
            save([FileBase '.pframp-' ElLoc '.tmp'],'out');%,'ThFrSpec','ThAmpSpec');
        end

        %         else %if file exists
        %             fprintf('File is computed, not overwriting\n');
        % %            load([FileBase '.pframp-' ElLoc],'-MAT');
        %         end
        out = CatStruct(out);
    case 'display'
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %display
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        load([FileBase '.pframp-' ElLoc],'-MAT');
        if ~exist('out','var')
            out = OutArgs;
        end
        [nSt, nClu] = size(out);

        for cnt=1:nClu
            figure(876)
            clf
            anyc=0;
            for w = 1:nSt %loop through states
                if isempty(out(cnt).El);  continue;    end
                anyc =anyc+1;
                fprintf('State: %s Cell: El=%d, Clu=%d, n=%d, log(Z)=%2.2f pR=%2.4f, R=%2.2f, k=%2.2f, pNU=%2.2f , V=%2.2f\n',...
                    States{w}, out(cnt).El, out(cnt).Clu, out(cnt).n, out(cnt).logZ, out(cnt).pR, out(cnt).R, out(cnt).k, out(cnt).pNU, out(cnt).V);


                %plot phase hist
                subplot(nSt,4,1+(w-1)*4)
                bar([out(cnt).phbin out(cnt).phbin+360],[out(cnt).phhist out(cnt).phhist]);
                axis tight
                ylabel(States{w});

                %plot shift modulation
                %subplot(nSt,4,2+(w-1)*4)

                %             plot(out(cnt).sh_lag,out(cnt).sh_logZ);
                %             axis tight
                %             yl = ylim; ylim([-1 yl(2)]);

                %             %plot shift params

                %             ax = plotyy(out(cnt).sh_frlag/eFs*1000,out(cnt).sh_frav,out(cnt).sh_amplag/eFs*1000,out(cnt).sh_ampav);
                %             hold(ax(1),'on'); hold(ax(2),'on');
                %             %plot(ax(1), out(cnt).sh_frlag,out(cnt).sh_frav+out(cnt).sh_frstd,'--');
                %             %plot(ax(2), out(cnt).sh_amplag,out(cnt).sh_ampav+out(cnt).sh_ampstd,'--');
                %             axis(ax(1),'tight');
                %             axis(ax(2),'tight');
                if 0
                    subplot(nSt*2,4,2+(w-1)*8)
                    smimage(out(cnt).t,linspace(out(cnt).fr_yb(1),out(cnt).fr_yb(end-1),10), out(cnt).ccg_fr,[smx 1],1,1);axis xy
                    xlim([-300 300]);
                    xlabel('CCG by freq');

                    subplot(nSt*2,4,6+(w-1)*8)
                    smimage(out(cnt).t,linspace(out(cnt).amp_yb(1),out(cnt).amp_yb(end-1),10), out(cnt).ccg_amp,[smx 1],1,1);axis xy
                    xlim([-300 300]);
                    xlabel('CCG by amp');
                end

                %plot fr. hist
                subplot(nSt,4,3+(w-1)*4)
                %        imagesc(out(cnt).fr_xb(1:end-1), out(cnt).fr_yb(1:end-1),unity(out(cnt).fr_hist)');
                smimage(out(cnt).fr_xb(1:end-1),out(cnt).fr_yb(1:end-1),out(cnt).fr_hist,[smx 1],1,1);axis xy
                xlabel('phase by freq.');
                axis tight; axis xy

                %plot time
                subplot(nSt,4,4+(w-1)*4)
                %imagesc(out(cnt).amp_xb(1:end-1), out(cnt).amp_yb(1:end-1),unity(out(cnt).amp_hist)');
                smimage(out(cnt).amp_xb(1:end-1),out(cnt).amp_yb(1:end-1),out(cnt).amp_hist,[smx 1],1,1);axis xy
                xlabel('phase by amp.');

            end
            suptitle([FileBase ' - ' num2str([out(cnt).El, out(cnt).Clu])]);
            drawnow
            if anyc>0%out(cnt).logZ>1.1 || out(cnt).V>1.7 % <0.0001 %| out(cnt).k>0.3

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
                    if curwin>5; curwin = curwin-5; end

                    switch sel
                        case 'normal'
                            break
                        case 'alt'
                            for curst=1:nSt
                                if isempty(out(c,curst).cID);  continue;    end
                                subplot(nSt,4,1+(curst-1)*4)
                                if curwin==2
                                    frax= linspace(out(c,curst).fr_yb(1),out(c,curst).fr_yb(end-1),10);
                                    [dummy curfr] = min(abs(frax-ymouse));
                                    bar(out(c,curst).t, out(c,curst).ccg_fr(:,curfr));
                                    axis tight
                                    xlim([-150 150]);

                                elseif curwin==3
                                    %[dummy curlag] = min(abs(out(c,curst).sh_lag-xmouse));
                                    %bar([out(c,curst).phbin out(c,curst).phbin+360],[out(c,curst).sh_hist(:,curlag); out(c,curst).sh_hist(:,curlag)]);
                                    ampax= linspace(out(c,curst).amp_yb(1),out(c,curst).amp_yb(end-1),10);
                                    [dummy curamp] = min(abs(ampax -ymouse));
                                    bar(out(cnt).t, out(c,curst).ccg_amp(:,curamp));
                                    axis tight
                                    xlim([-150 150]);

                                elseif curwin==4
                                    [dummy curfr] = min(abs(out(c,curst).fr_yb(1:end-1)-ymouse));
                                    bar(out(c,curst).fr_xb(1:end-1), out(c,curst).fr_hist(:,curfr));
                                    axis tight

                                elseif curwin==5
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
                        out(cnt).cID, out(cnt).logZ, out(cnt).pR, out(cnt).R, out(cnt).k, out(cnt).pNU, out(cnt).V);

                    reportfig(gcf, 'PhaseAnal',0,[FileBase '  ' rstr ],150);
                end
            end

        end
    case 'display2'
        if ~FileExists([FileBase '.' mfilename '.mat'])
            load([FileBase '.pframp-' ElLoc],'-MAT');
        else
            load([FileBase '.' mfilename '.mat']);
        end
        if ~exist('out','var')
            out = OutArgs;
        end
        CluLoc = load([FileBase '.cluloc']);
        [nSt, nClu] = size(out);
        figure(3232);clf
        for cnt=1:nClu
            if isempty(out(cnt).El) continue; end
            subplotfit(cnt,nClu);
            bar([out(cnt).phbin out(cnt).phbin+360],[out(cnt).phhist out(cnt).phhist]);
            axis tight
            mloc = CluLoc(find(CluLoc(:,1)==out(cnt).El & CluLoc(:,2)==out(cnt).Clu),3);
            title(num2str([out(cnt).El out(cnt).Clu mloc]));
              fprintf('State: %s Cell: El=%d, Clu=%d, n=%d, log(Z)=%2.2f pR=%2.4f, R=%2.2f, k=%2.2f, pNU=%2.2f , V=%2.2f\n',...
                    out(cnt).State, out(cnt).El, out(cnt).Clu, out(cnt).n, out(cnt).logZ, out(cnt).pR, out(cnt).R, out(cnt).k, out(cnt).pNU, out(cnt).V);

        end

    case 'group'
        if FileExists([FileBase '.pframp-' ElLoc])
            load([FileBase '.pframp-' ElLoc],'-MAT');
            %out = OutArgs;
        else
            load([FileBase '.' mfilename '.mat']);
        end
            if ~isempty(out)
                out= CatStruct(out,{'FileBase','El','Clu',...
                    'State','pR','R','k','logZ','th0','n','StateTime','V','pNU'});
            end


end

        %Out = CatStruct(out,[],3);





