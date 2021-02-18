%function out = TrigPtSpec(FileBase,Trig, ElLoc,States, Overwrite, Display, Recompute, LfpChan,Win,Step,Duration, Subset)
function out = TrigPtSpec(FileBase,varargin)
[Trig, ElLoc, States,         Overwrite,  Display, Recompute, LfpChan,  Win,    Step,  Duration, Subset] = DefaultArgs(varargin, ...
    {[], 'c',   {'REM','RUN'},  0,          0,       0,         [],     50,     10,    1,      1           });

out = struct([]);
Par = LoadPar([FileBase '.xml']);
MinPeriod = 10; %seconds for theta periods selection
if isstr(ElLoc(1))
    El = find(strcmp(Par.ElecLoc,ElLoc));
else
    El = ElLoc;
    ElLoc = num2str(El);
end
if ischar(States); States = {States}; end;
nSt = length(States);
FreqRange =[5 12];
eFs = 1250;
frThr = 7; %Hz
MinSpikes =2; %hoping to improve the specgrams at high freq.


if (~FileExists([FileBase '.TrigPtSpec_' ElLoc '.mat']) | Overwrite) & Recompute

    [Res,Clu,nclu,sph,elclu] = ReadEl4CCG(FileBase,El);
    if isempty(Res)
        fprintf('no units in %s in this file \n',ElLoc);
        return;
    end

    Res = round(Res/16)+1;

    %     if isstr(Trig)
    %         if ~isempty(strfind(Trig,'evt'))
    %         Trig =
    if iscell(Trig)
        [EvtRes, EvtClu, EvtLabels, Labels] = LoadEvt([FileBase '.' Trig{1} '.evt'],eFs, Trig{2});
        Trig = EvtRes;
    end
    out = struct([]);
    cnt =1;
    for w = 1:length(States) %loop through states

        %load theta periods (more then MinPeriod)
        %        if FileExists([FileBase '.sts.' States{w}])
        Period = SelectStates(FileBase , States{w}, MinPeriod*eFs,0);
        %       else
        %          fprintf('No %s periods. Return empty output\n',States{w});
        %          continue;
        %      end
        if isempty(Period);
            fprintf('No %s periods. Return empty output\n',States{w});
            continue;
        end
        PeriodTime = sum(diff(Period,[],2))/eFs; % in seconds
        fprintf('\nFile %s has %2.1f seconds of %s \n',FileBase,PeriodTime,States{w});

        [InRes In]= SelectPeriods(Res,Period,'d',1);
        InClu = Clu(In);
        uInClu= unique(InClu);
        InTrig = SelectPeriods(Trig,Period,'d',1);

        FirRate = FiringRate(InRes,InClu,uInClu,eFs);
        GoodCells = uInClu(find(FirRate > frThr & ~isnan(FirRate)));%clusters of cells that fire enough

        goodRes = InRes(ismember(InClu,GoodCells));
        goodClu = InClu(ismember(InClu,GoodCells));
        myEl = elclu(uInClu(GoodCells),1);
        myClu = elclu(uInClu(GoodCells),2);
        out(cnt).El = myEl;
        out(cnt).Clu = myClu;
        out(cnt).FirRate = FirRate(GoodCells);

        fprintf('Processing cells %d %d\n',myEl, myClu);
        if isempty(LfpChan)
            CluLoc = ClusterLocation(FileBase, myEl);
            LfpCh = CluLoc(CluLoc(:,1)==myEl(1) & CluLoc(:,2)==myClu(1),3);
        end

        eeg = LoadBinary([FileBase '.eeg'],LfpChan,Par.nChannels);
        eeg = WhitenSignal(eeg);
        %eeg = eeg(:);
        step = 2^nextpow2(eFs*Step/1000);
        win=2^nextpow2(eFs*Win/1000);
        %            for cnt=1:400
        [out(cnt).y, out(cnt).f, out(cnt).t, out(cnt).phi, out(cnt).yerr, out(cnt).phierr, out(cnt).ynorm, out(cnt).ycom, out(cnt).phinorm ]...
            = mtptchgram_trig(struct('Eeg',eeg,'Res',goodRes,'Clu',goodClu),...
            InTrig(1:Subset:end), Duration*eFs , struct('nFTFT',2^10,'WinLength',win, 'WinStep',step,'NW',1,'MinSpikes',MinSpikes,'FreqRange',[20 100]));
        %           end
        [out(cnt).ccg,out(cnt).tbin] = Trains2CCG({InTrig,goodRes},{1,goodClu},Step,round(win/Step),1250);
        out(cnt).State = States{w};
        cnt = cnt+1;


    end
    if    ~FileExists([FileBase '.TrigPtSpec_' ElLoc '.mat']) | Overwrite
        save([FileBase '.TrigPtSpec_' ElLoc '.mat'],'out');
    end

end


if Display

    if ~Recompute
        if FileExists([FileBase '.TrigPtSpec_' ElLoc '.mat'])
            load([FileBase '.TrigPtSpec_' ElLoc '.mat']);
        elseif FileExists([FileBase '.TrigPtSpec.mat'])
            load([FileBase '.TrigPtSpec.mat']);
            out = OutArgs;
        else
            error('have nothing to display');
        end
    end
    mode = 1;
    switch mode
        case 1
            t = out(1).t;    f = out(1).f; myf = find(f>20 & f<100);
            nch = size(out.y,3);
            for i=1:nch
                subplot(nch,nch,(i-1)*nch+i);
                npow = sq(out(1).y(:,myf,i,i)).*repmat(f(myf).^2',length(t),1);
                imagesc(t*1000,f(myf),npow');axis xy
                for j=i+1:nch
                    subplot(nch,nch,(i-1)*nch+j);
                    imagesc(t*1000,f(myf),sq(out(1).y(:,myf,i,j))');axis xy
                    title(['coherence ' num2str(i) ' to ' num2str(j)]);
                end
            end

        case 2

            nCells = length(out);
            t = out(1).t;    f = out(1).f; myf = find(f>30 & f<100);
            ni=4; nj=2;
            for i=1:nCells

                figure(2126);
                clf
                subplot(ni,nj, 1)
                %        npow = log(sq(out(i).y(:,myf,1,1)));
                npow = sq(out(i).y(:,myf,1,1)).*repmat(f(myf).^2',length(t),1);
                imagesc(t*1000,f(myf),npow');axis xy
                title('power of eeg');

                subplot(ni,nj, 2)
                npowspk = sq(out(i).y(:,myf,2,2)) ./repmat(interp1(out(i).tbin,sq(out(i).ccg(:,1,2)), out(i).t),length(myf),1)';
                imagesc(t*1000,f(myf),npowspk');axis xy
                title('power of spikes');

                subplot(ni,nj, 3)
                imagesc(t*1000,f(myf),atanh(sq(out(i).y(:,myf,1,2)))');axis xy
                title('coherence spks to eeg');
                hold on
                phvar = 1./(sq(out(i).phierr(:,myf(1:5:end),1,2))+eps);
                u = real(phvar.*exp(sqrt(-1)*sq(out(i).phi(:,myf(1:5:end),1,2))));
                v = imag(phvar.*exp(sqrt(-1)*sq(out(i).phi(:,myf(1:5:end),1,2))));
                quiver(t*1000,f(myf(1:5:end)),u',v',1, 'Color','k');

                subplot(ni,nj, 4)
                bar(out(i).tbin, out(i).ccg(:,1,2));axis tight
                yl = ylim;    ylim([mean(yl) 1.3*yl(2)]);
                title('spike ccg');

                %plot again with normalization
                subplot(ni,nj, 5)
                imagesc(t*1000,f(myf),atanh(sq(out(i).ynorm(:,myf,1,2)))');axis xy
                title('norm coherence spks to eeg');
                hold on
                phvar = 1./(sq(out(i).phierr(:,myf(1:5:end),1,2))+eps);
                u = real(phvar.*exp(sqrt(-1)*sq(out(i).phinorm(:,myf(1:5:end),1,2))));
                v = imag(phvar.*exp(sqrt(-1)*sq(out(i).phinorm(:,myf(1:5:end),1,2))));
                quiver(t*1000,f(myf(1:5:end)),u',v',1, 'Color','k');


                subplot(ni,nj, 6)
                imagesc(t*1000,f(myf),abs(sq(out(i).ycom(:,myf,1,2)))');axis xy
                title('average coherence spks to eeg');
                hold on

                subplot(ni,nj, 7)
                imagesc(t*1000,f(myf),(sq(out(i).phinorm(:,myf,1,2)))');axis xy
                title('norm phase of spks to eeg');

                subplot(ni,nj, 8)
                imagesc(t*1000,f(myf),(sq(out(i).phi(:,myf,1,2)))');axis xy
                title('phase of spks to eeg');


                tit = [FileBase ' : El=' num2str(out(i).El) ', Clu=' num2str(out(i).Clu) ...
                    ' (' out(i).State '), Rate=' num2str(out(i).FirRate)];
                suptitle(tit);


                if 1
                    [x,y,b] = PointInput(1);
                    if b>2
                        reportfig(gcf,'TrigGammaSpec',0,FileBase,150);
                    end
                end

            end
    end


end
