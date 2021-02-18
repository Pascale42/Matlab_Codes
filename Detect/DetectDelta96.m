%function out = DetectDelta96(FileBase, fMode, Channels,FreqRange)
% fMode = 'channel', 'compute', 'sort'
function out = DetectDelta96(FileBase, varargin)
Par = LoadPar([FileBase '.xml']);
eFs= 1250;

[fMode, Channels, FreqRange] = DefaultArgs(varargin,{'detect',load([FileBase '.ecsink']),[100 600]});
State ='SWS';
%if isempty(RepChannels)
[RepChan Tetrode] = RepresentChan(FileBase);
%end

nT = FileLength([FileBase '.eeg'])/Par.nChannels/2;
%parameters set:

NegPosLag = round(eFs*0.3);
ResampleTo = 25; %Hz
swin = 50;
rs=10; % resamplng factor
reFs = eFs/rs;
Period = load([FileBase '.sts.' State]);


switch fMode
    case 'compute'

        eeg = LoadBinary([FileBase '.csd'], Channels, Par.nChannels);
        %csd = LoadBinary([FileBase '.csd'], Channels, Par.nChannels, [], [], 'double', Period, 1)';

        slowfeeg = ButFilter(eeg,2,[2 7]/(eFs/2),'bandpass');
        weeg = WhitenSignal(eeg);
        fastfeeg = ButFilter(weeg,2,FreqRange/(eFs/2),'bandpass');
        gpow = Filter0(ones(swin,1)/swin,abs(fastfeeg));
        gpow = ButFilter(abs(fastfeeg),2,20/625,'low');

        
        [Chans CxEl] = GetChannels(FileBase,'c');
        [Res,Clu,nclu] = ReadEl4CCG(FileBase,CxEl);
        Res = round(Res/16);
        [Res ind] =SelectPeriods(Res,Period,'d',1);
        Clu =Clu(ind);
        [uClu dummy Clu] = unique(Clu);
        nClu = max(Clu);

        %detect onset of upstates
        [Burst, BurstLen, SpkPos, OutOf] = SplitIntoBursts(Res, 1.25*150);
        
        [out.trw.tav out.trw.tstd out.trw.twav out.trw.f out.trw.t] = TrigScalogram(eeg,Res(Burst),0.5,1250,[1 650]);
        gseg = GetSegs(gpow,Res(Burst)-625,1250,[]);
        [u s v]= svd(gseg,0);
        out.gsvd.s = diag(s);
        out.gsvd.u = u(:,1:5);
        out.gsvd.v = v(:,1:5);

    case 'detect'

        eeg = LoadBinary([FileBase '.csd'], Channels, Par.nChannels);
        %csd = LoadBinary([FileBase '.csd'], Channels, Par.nChannels, [], [], 'double', Period, 1)';

        slowfeeg = ButFilter(eeg,4,10/(eFs/2),'low');
        weeg = WhitenSignal(eeg);
        fastfeeg = ButFilter(weeg,4,FreqRange/(eFs/2),'bandpass');

        gpow = Filter0(ones(swin,1)/swin,abs(fastfeeg));
        gpow = smooth(gpow,50,'lowess');
        %        gpow = ButFilter(gpow,4,10/(eFs/2),'low');
        
        if 0
        
        [Chans CxEl] = GetChannels(FileBase,'c');
        [Res,Clu,nclu,dummy,ElByCluCx] = ReadEl4CCG(FileBase,CxEl);
        Res = round(Res/16);
        [Res ind] =SelectPeriods(Res,Period,'d',1);
        Clu =Clu(ind);
        [uClu dummy Clu] = unique(Clu);
        nClu = max(Clu);
        %detect onset of upstates
        [Burst, BurstLen, SpkPos, OutOf] = SplitIntoBursts(Res, 1.25*150);
        end
        %         rgpow = resample(gpow,1,rs);
        %         rfeeg = resample(slowfeeg,1,rs);
        %         rPeriod = round(Period/rs);

        len = round(100*eFs/1000);
        win = [-ones(len,1); ones(len,1)];
        wins = win.*hamming(len*2);
        wins=wins./length(wins);

        dgam = Filter0(wins, gpow);

        gbeg = LocalMinima(-dgam,eFs*200/1000,0);
        gend = LocalMinima(dgam,eFs*200/1000,0);

        MinISI = round(eFs*0.2);
        thr = 0.5*std(SelectPeriods(slowfeeg,Period,'c',1));
        pkspos = LocalMinima(-slowfeeg, MinISI,-thr);
        pkspos = SelectPeriods(pkspos,Period,'d',1);

        pksneg = LocalMinima(slowfeeg);
        pksneg = SelectPeriods(pksneg,Period,'d',1);

        [negr pksnir pkspi] = NearestNeighbour(pksneg, pkspos,'right',150*eFs/1000);
        pkspos = pkspos(pkspi);
        [negl pksnil pkspi] = NearestNeighbour(pksneg, pkspos,'left',250*eFs/1000);
        pkspos = pkspos(pkspi);
        negr = negr(pkspi);

        gammin = LocalMinima(gpow,10*eFs/1000);
        [gammin gmi pkspi] = NearestNeighbour(gammin,pkspos,'both',50*eFs/1000);
        pkspos = pkspos(pkspi);
        negr = negr(pkspi);
        negl = negl(pkspi);

        gammax = LocalMinima(-gpow,10*eFs/1000);
        [gammax gmi pkspi] = NearestNeighbour(gammax,negr,'both',50*eFs/1000);

        n2ntau = negr-negl;
        posamp = slowfeeg(pkspos);
        nlamp = slowfeeg(negl);
        nramp = slowfeeg(negr);
        n2pamp = min(posamp-nlamp,posamp-nramp);
        posd2 = posamp./n2ntau;
        pos2gam = pkspos-gammin;
        gamnr = gpow(negr);
        gamnl = gpow(negl);
        gampos = gpow(gammin);
        n2pgam = gamnr-gampos;
        
        
%        keyboard

        %now find Bursts close to trough
%        [NearBursts nbi] = NearestNeighbour(negr,Res(Burst),'both',30*eFs/1000);

        twin = round(300*eFs/1000/10);
        segs = GetSegs(resample([eeg gpow],1,10),round(negr/10)-twin,2*twin+1,[]);
        nsegs = segs./repmat(std(segs,[],2),[1 size(segs,2) 1]);

%         sgpow = ButFilter(gpow,4,20/625,'low');
%         
%         gmin = LocalMinima(sgpow);
%         gmax = LocalMinima(-sgpow);
%         gmin = SelectPeriods(gmin,Period,'d',1);
%         gmax = SelectPeriods(gmax,Period,'d',1);
%         [gmaxr dummy gmini] = NearestNeighbour(gmax,gmin,'right');
%         gmin= gmin(gmini);
%         [gmaxl dummy gmini] = NearestNeighbour(gmax,gmin,'left');
%         gmin= gmin(gmini);
%         gmaxr = gmaxr(gmini);
%         gminlen = gmaxr-gmaxl;
%         gi = gminlen<500;
%         gminlen =gminlen(gi);
%         gminamp = [sgpow(gmaxl(gi)) sgpow(gmaxr(gi))];
%         
         twin = round(50*eFs/1000);
         gsegs = GetSegs(gpow,round(pkspos)-twin,2*twin+1,[]);
         mgampos = median(gsegs);
         vgampos = std(gsegs);
%         keyboard
        
        fet = [n2ntau posamp -nramp n2pamp posd2 pos2gam -gampos gamnr n2pgam -mgampos' -vgampos'];
        labels= {'n2ntau' 'posamp' 'nramp' 'n2pamp' 'posd2' 'pos2gam' 'gampos' 'gamnr' 'n2pgam' 'mgampos' 'vgampos'};
%        ExamSegs('init', nsegs, fet, labels, [FileBase '.delta.clu']);
        
        Kluster([FileBase '.delw'],{repmat(fet,[1 1 2]),permute(nsegs,[3 1 2]),negr},[],1250);
        return       
%       keyboard
        
        out.fet = fet;
        out.labels = labels;
        out.negr = negr;
        out.pkspos = pkspos;
        out.gammax = gammax;
        
%        segs = GetSegs([eeg gpow],round(negr)-650,1251,[]);

        % %         [u s v] = svd(segs,0);
        % %         fet = [posamp negamp gdepth rgpow(pksneg) pksneg-pkspos v(:,1:3)];
        % %         labels = {'posamp','negamp', 'gamdepth','posgam','tau-pos-neg','pc1','pc2','pc3'};
        %
        figure(2121);clf
        for i=1:9
            subplot(3,3,i)
            bs=linspace(min(fet(:,i)),max(fet(:,i)),20);
            h1=hist(fet(:,i),bs);
            h2=hist(fet(nbi,i),bs);
            plot(bs,unity([h1(:) h2(:)]));axis tight
            title(labels{i});
        end
        reportfig(gcf,'detect_delta',0,[FileBase ' - calibration'],200,0);
        return



        [ledge leind  ] = NearestNeighbour(gbeg, pkspos);
        [redge reind  ] = NearestNeighbour(gend, pkspos);


        ta = [1:length(rfeeg)]/1250;
        %        PlotTraces(unity([rfeeg rgpow srgpow dgam]),ta,reFs,3);
        PlotTraces(unity([rfeeg rgpow srgpow]),ta,reFs,3);
        hold on
        Lines(ta(pkspos),[],'r');
        %         Lines(ta(redge),[],'b');
        %         Lines(ta(ledge),[],'g');
        Lines(Res(Burst)/10/1250,[],'k');
        TimeBrowse(100,100)
        keyboard

        gamsegs = GetSegs(unity(rgpow),pkspos-30,60,[]);
        pksgam = mean(gamsegs);


        %     thr = 1*std(gpow);
        %     gamtrough = LocalMinima(gpow, MinISI, -thr);

        %now find gamma transitions around delta peaks


        res = [Delta.neg.t(:); Delta.pos.t(:)];

        [res, ind ] = sort(res);
        lab = lab(ind);
        MakeEvtFile(res,[FileBase '.del.evt'],{'1','2'},eFs);

    case 'cluster'
    load([FileBase '.' mfilenname '.mat']);
    

    case 'display'
        load([FileBase '.' mfilename '.mat']);

        t = OutArgs.trw.t; t=t-median(t);
        f = OutArgs.trw.f;
        tav = OutArgs.trw.tav;
        ti = find(abs(t)<0.4);
        figure(21211);clf
        subplot(221)
        imagelog(t(ti),f,log(tav(ti,:)'));

        subplot(222)
        imagelog(t(ti),f,tav(ti,:)'.*repmat(f'.^2,length(ti),1)');

        subplot(223)
        imagelog(t(ti),f,log(tav(ti,:)')-log(repmat(mean(tav(ti,:)),length(ti),1)'));

        subplot(224);cla
        gev = sign(median(OutArgs.gsvd.v(:,1)))*OutArgs.gsvd.u(:,1);
        plot([-625:624]/1.25,gev);axis tight
        hold on
        Lines(0,[],'r');

        % keyboard
        return;



end