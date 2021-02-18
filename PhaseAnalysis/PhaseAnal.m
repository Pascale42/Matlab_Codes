function out = PhaseAnal(FileBase,varargin)
%function out = PhaseAnal(FileBase,State,ElLoc,Display,Overwrite, Report,UseDent)
%does all the phase estimation for all cells from certian location
% UseDent - flag to use dentate phase if it exists
[States, ElLoc,Display, Overwrite, Report] = DefaultArgs(varargin,{{'REM','RUN'},'c',1,0,0});

MinPeriod = 1; %seconds for theta periods selection
MinNumSpks = 100;
ThCh=1;
par = LoadPar([FileBase '.xml']);
sFs = 1e6/par.SampleTime;
eFs = 1250;%GetEegFs(FileBase);
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
if ( ~FileExists([FileBase '.pha']) | Overwrite) & ~Display

    %channels used for phase calculation
    Channels = load([FileBase '.eegseg.par']);
    Channels = Channels(1);
    nCh = length(Channels);

    %now load all spikes
    [Res, Clu, nClu, spikeph, ClustByEl, cID] = ReadEl4CCG(FileBase,El);
    %keyboard
    clear spikeph; %don't use the old one
    Res =  max(1,round(Res*eFs/sFs)); %get spike times to the sample rate of phase.avoid zero

    %now load the phases vectors (and amplitude)
    [ThPh, ThAmp, ThFr] = ThetaParams(FileBase);
    ThPhu = MakeUniformDistr(ThPh(:,ThCh));
    
    out = struct([]);
    cnt=1;
    for where = 1:length(States) %loop through states

        %load theta periods (more then MinPeriod)
        if FileExists([FileBase '.sts.' States{where}])
            AllPeriod = load([FileBase '.sts.' States{where}]);
            Period = AllPeriod(find(diff(AllPeriod,1,2)>eFs*MinPeriod),:);
            Period = Period*eFs;
        else
            fprintf('No %s periods. Return empty output\n',States{where});
            continue;
        end
        if isempty(Period);
            fprintf('No %s periods. Return empty output\n',States{where});
            continue;
        end

        PeriodTime = sum(diff(Period,[],2))/eFs; % in seconds
        fprintf('\nFile %s has %2.1f seconds of %s \n',FileBase,PeriodTime,States{where});

        Focus = WithinRanges(Res, Period);
        Focus = Focus==1; % make it logical to use fast subscription

        %select spikes that are in that State
        %         uClu = unique(Clu(Focus));
        %         nClu = length(uClu);
        cnt = 1;
        %         myClu = ismember([1:nClu],uClu); %indexes of clusters that come into cohernce computation (present in this state)
        fprintf('State : %s - PROCESSING\n',States{where});


        % fields of struct PhOut :  pR, th0, R, phconfint, mu, k, pNU, V, phhist, phbin
        %now loop through cells
        for c=1:nClu
         
            %limit Res to the Periods
            myRes = Res(find(Clu==c & Focus));
            if length(myRes)<MinNumSpks
                continue;
            end
            out(c,where).cID = cID(c);
            out(c,where).ElClu = ClustByEl(c,:);
            myph = ThPh(myRes,ThCh);
            myphu = ThPhu(myRes);
            myfr = ThFr(myRes,ThCh);
            myamp = log(ThAmp(myRes,ThCh));
            % test for unimodal vs uniform distr
%            [out(c,where).pR out(c,where).th0 out(c,where).R out(c,where).logZ] = RayleighTest(myph);
            [out(c,where).RT] = RayleighTest(myph);
            
            % test for any alternative to unofrm
            [H, out(c,where).V, out(c,where).pNU] = TestPhaseNonuniform(myph);

            %out(c,where).logZ = log(-log(out(c,where).pval));
            out(c,where).th0 = out(c,where).th0*180/pi;

            [cm, r] = CircConf(myph,0.01,200);
            out(c,where).phconfint = r*180/pi;
            [out(c,where).mu out(c,where).k] = VonMisesFit(myph);

            out(c,where).n = length(myRes);

            %cvPh = PhaseFit(myph);
            %out(c,where).cvL1 = cvPh.L1;
            %out(c,where).cvL2 = cvPh.L2;
            myph = myph*180/pi;
            %myphall = mod(myphall*180/pi,360);

            edges = [-180:20:180];
            [out(c,where).phhist ] = histcI(myph,edges);
            out(c,where).phbin = (edges(2:end)+edges(1:end-1))/2;
            out(c,where).phhist = out(c,where).phhist';

            %compute lagged phase modulation

            MaxShift = round(pFs*2000/1000); BinSize = round(pFs/1000*20);
            [out(c,where).sh_p, out(c,where).sh_lag, out(c,where).sh_th0, out(c,where).sh_r, out(c,where).sh_logZ, out(c,where).sh_hist] = ...
                ShiftPhase(ThPh(:,ThCh),myRes,MaxShift,BinSize,Period);
            out(c,where).sh_lag =out(c,where).sh_lag/pFs*1000;

            %compute shifted freq.
            [out(c,where).sh_frav,out(c,where).sh_frlag, out(c,where).sh_frstd] = ShiftParam(ThFr(:,ThCh),myRes,MaxShift,BinSize,Period);
            % .... and amplitude of theta
            [out(c,where).sh_ampav,out(c,where).sh_amplag, out(c,where).sh_ampstd] = ShiftParam(ThAmp(:,ThCh),myRes,MaxShift,BinSize,Period);


            %compute phase histogram as a function of freq.
            prct = prctile(myfr,[5 95]); gi = find(myfr>prct(1) & myfr<prct(2));
            [out(c,where).fr_hist, out(c,where).fr_xb,out(c,where).fr_yb] = hist2([myph(gi) myfr(gi); myph(gi)+360 myfr(gi)], 18, 10);

            %compute phase histogram as a function of freq.
            prct = prctile(myamp,[5 95]); gi = find(myamp>prct(1) & myamp<prct(2));
            [out(c,where).amp_hist, out(c,where).amp_xb,out(c,where).amp_yb] = hist2([myph(gi) myamp(gi); myph(gi)+360 myamp(gi)], 18, 10);


            fprintf('State: %s Cell %10f : n=%d, log(Z)=%2.2f pR=%2.4f, R=%2.2f, k=%2.2f, pNU=%2.2f , V=%2.2f\n',...
                States{where}, out(c,where).cID, out(c,where).n, out(c,where).logZ, ...
                out(c,where).pR, out(c,where).R, out(c,where).k, out(c,where).pNU, out(c,where).V);

        end

        %Out(where) = a;
        %a=CatStruct(out);

    end
    save([FileBase '.pha'],'out');
else
    load([FileBase '.pha'],'-MAT');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%display
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if Display
    [nClu, nSt] = size(out);
    mode =1;
    switch mode
        case 1
            for where = 1:nSt %loop through states
                 if isempty(cat(1,out(:,where).cID))
                     continue;
                 end
                figure(876)
                clf
                for c=1:nClu
                    if isempty(out(c,where).cID);  continue;    end
                    
                    subplotfit(c,nClu);
                    bar(out(c,where).phbin,out(c,where).phhist);hold on
                    bar(out(c,where).phbin+360,out(c,where).phhist);
                    axis tight
                    title([FileBase '  ' num2str(out(c,where).ElClu(1)) '/' num2str(out(c,where).ElClu(2)) ' - ' States{where} ' : ' num2str(out(c,where).logZ)]);
                   %suptitle(FileBase);
                    
                end
                hto = round(rand(1,1)*1000)+1; %copy to axis
%                prtstring = [FileBase '_' num2str(out(c,where).ElClu(1)) '_' num2str(out(c,where).ElClu(2)) '.eps'];
                prtstring = [FileBase '_' States{where} '_' num2str(hto)];
                copyax('on',hto,gcf,prtstring);
               %                 figure(hto);
%                 prteps();
%                 %waitforbuttonpress;
%                 figure(876);
                pause
                copyax('off');
            end
        case 2


            for c=1:nClu
                figure(876)
                clf


                anyc=0;
                for where = 1:nSt %loop through states
                    if isempty(out(c,where).cID);  continue;    end


                    anyc =anyc+1;
                    %if out(c,where).logZ<1.5 continue; end
                    fprintf('State: %s Cell: El=%d, Clu=%d, n=%d, log(Z)=%2.2f pR=%2.4f, R=%2.2f, k=%2.2f, pNU=%2.2f , V=%2.2f\n',...
                        States{where}, out(c,where).ElClu(1), out(c,where).ElClu(2), out(c,where).n, out(c,where).logZ, out(c,where).pR, out(c,where).R, out(c,where).k, out(c,where).pNU, out(c,where).V);


                    %plot phase hist
                    subplot(nSt,4,1+(where-1)*4)
                    bar([out(c,where).phbin out(c,where).phbin+360],[out(c,where).phhist out(c,where).phhist]);
                    axis tight
                    ylabel(States{where});

                    %plot shift modulation
                    %subplot(nSt,4,2+(where-1)*4)
                    subplot(nSt*2,4,2+(where-1)*8)
                    plot(out(c,where).sh_lag,out(c,where).sh_logZ);
                    axis tight
                    yl = ylim; ylim([-1 yl(2)]);

                    %plot shift params
                    subplot(nSt*2,4,6+(where-1)*8)
                    ax = plotyy(out(c,where).sh_frlag/pFs*1000,out(c,where).sh_frav,out(c,where).sh_amplag/pFs*1000,out(c,where).sh_ampav);
                    hold(ax(1),'on'); hold(ax(2),'on');
                    %plot(ax(1), out(c,where).sh_frlag,out(c,where).sh_frav+out(c,where).sh_frstd,'--');
                    %plot(ax(2), out(c,where).sh_amplag,out(c,where).sh_ampav+out(c,where).sh_ampstd,'--');
                    axis(ax(1),'tight');
                    axis(ax(2),'tight');

                    %plot fr. hist
                    subplot(nSt,4,3+(where-1)*4)
                    imagesc(out(c,where).fr_xb(1:end-1), out(c,where).fr_yb(1:end-1),out(c,where).fr_hist');
                    axis tight; axis xy

                    %plot time
                    subplot(nSt,4,4+(where-1)*4)
                    imagesc(out(c,where).amp_xb(1:end-1), out(c,where).amp_yb(1:end-1),out(c,where).amp_hist');
                    axis tight; axis xy

                end
                drawnow
                if anyc>0%out(c,where).logZ>1.1 || out(c,where).V>1.7 % <0.0001 %| out(c,where).k>0.3

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
                            case 'aaaalt'
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
                            case 'alt'
                                reportfig(gcf, 'PhaseAnal',0,[FileBase '  ' rstr ],150);

                                %keyboard
                        end
                    end

                    if Report
                        rstr = sprintf('Cell %10f : log(Z)=%2.2f pR=%2.4f, R=%2.2f, k=%2.2f, pNU=%2.2f , V=%2.2f\n',...
                            out(c,where).cID, out(c,where).logZ, out(c,where).pR, out(c,where).R, out(c,where).k, out(c,where).pNU, out(c,where).V);

                        reportfig(gcf, 'PhaseAnal',0,[FileBase '  ' rstr ],150);
                    end


                end %of for where

            end%of cells loop
    end %of switch
end

%Out = CatStruct(out,[],3);





