%function Power2Units(FileBase, Electrodes, Window, FreqRange, States)
function out = Power2Units(FileBase,varargin)
[fMode, Electrodes, Window, FreqRange,States] = DefaultArgs(varargin,{'compute','c',1,[1 100],{'REM','RUN'}});

Par = LoadPar([FileBase '.par']);
SampleRate = 1e6/Par.SampleTime;

if isstr(Electrodes)
    Electrodes = find(strcmp(Par.ElecLoc,Electrodes));
end

switch fMode
    case {'compute','slow','fast','adjtrig'}

        [Res, Clu, Map ] = LoadCluRes(FileBase, Electrodes);
        [Res, Clu] = RemoveBursts(Res,Clu,10*Par.SampleRate/1000);
        Channels = RepresentChan(FileBase);
        load([FileBase '.thpar.mat.'],'ThPh');
        wcnt=0;
        for w=1:length(States)
            if ~FileExists([FileBase '.sts.' States{w}])
                continue;
            else
                Period = load([FileBase '.sts.' States{w}]);
                wcnt=wcnt+1;
            end
            out(wcnt).State = States{w};
            
            [myRes ind] = SelectPeriods(round(Res*1250/SampleRate), Period,'d',1,1);
            myClu = Clu(ind);
            [uClu dummy myClu] = unique(myClu);
            nClu = length(uClu);
            out(wcnt).Map = Map(uClu,2:3);
            out(wcnt).nSpk = Accumulate(myClu,1);
            out(wcnt).Channels = Channels;
            
            if strcmp(fMode,'slow') | strcmp(fMode,'compute')
                Window = 2^round(log2(Window*25)); % now in eeg sample rate
                nFFT = Window*2; Fs =1250/25; nOverlap = Window - Window/4; NW= 2;
                %            MaxFreq =100;
                arMod = [];pow=[];
                for i=1:length(Channels)
                    Eeg = LoadBinary([FileBase '.eeg'], Channels(i), Par.nChannels,Period);
                    Eeg = resample(Eeg,1,25);
                    % Eeg = WhitenSignal(Eeg);
                    [Eeg, arMod] = WhitenSignal(Eeg,[],1,arMod);
                    [pow(:,:,i), f, t]=mtcsglong(Eeg, nFFT,Fs,Window,nOverlap,NW,'linear');
                end
                t = t+Window/2;
                [dummy, tRes , gspk] =NearestNeighbour(t,myRes/1250,'both');

                Firing = accumarray([tRes, myClu],1,[length(t) max(myClu)]);
                Firing=Firing/Window*50;%normalize firing rate to Hz
                mpow = mean(pow);
                for c=uClu'
                    if sum(myClu==c)>5
                        out(wcnt).Slow.AvPow(:,:,c) = log10(sq(mean(pow(tRes(myClu==c),:,:))))-sq(log10(mpow));%check

                        for fi = 1:length(f)
                            for chi =1:length(Channels)
                                [out(wcnt).Slow.Rho(fi,chi,c) out(wcnt).Slow.Sig(fi,chi,c)] = RankCorrelation(log10(pow(:,fi,chi)),Firing(:,c));%check
                            end
                        end
                    end
                end
                out(wcnt).Slow.f = f;
                if 0
                    pow = log(reshape(pow,size(pow,1),[]));
                    [LU, LR, FSr , VT] = erpPCA(pow);
                    nPCA=30;
                    out(wcnt).Slow.eigval = VT(1:nPCA,:);
                    out(wcnt).Slow.Score = FSr(:,1:nPCA);

                    out(wcnt).Slow.eigvec = reshape(LR(:,1:nPCA),[length(f), length(Channels), nPCA]);
                end
                % now gamma
                save([FileBase '.slowspec.' States{w}],'pow','f','t');
                clear Eeg y f t myRes;
            end
            
            if strcmp(fMode,'fast') | strcmp(fMode,'compute') 
                [Spec t f] = LoadSpecs(FileBase,States{w}); % t x f x ch
                Spec = exp(Spec);
                out(wcnt).Fast.f=f;
                mSpec = sq(mean(Spec));
                [myRes ind] = SelectPeriods(round(Res*1250/SampleRate), Period,'d',1);

                [dummy, tRes , gspk] =NearestNeighbour(t,myRes/1250,'both');

                phedge = linspace(-pi,pi,9);
                [dummy phind ] = histcI(ThPh(myRes),phedge);
                out(wcnt).Fast.phbin = (phedge(2:end)+phedge(1:end-1))*180/pi/2;
    
                out(wcnt).Fast.tbin = (t(2)-t(1))*[-200:200];
                % this is now: freq x ch x tlag x cell
                for c=uClu'
                    if sum(myClu==c)>2
                        for fi = 1:length(f)
                            for chi =1:length(Channels)
                                out(wcnt).Fast.TrAvPow(fi,chi,:,c) = log10(TriggeredAv(sq(Spec(:,fi,chi)),200,200,tRes(myClu==c)))...
                                    - log10(mSpec(fi,chi));
                            end
                        end
                        for j=1:8
                            if sum(myClu==c & phind==j)>5
                                out(wcnt).Fast.Pow2Ph(:,:,c,j) = log10(sq(mean(Spec(tRes(myClu==c & phind==j),:,:))))-log10(sq(mSpec));
                            end
                        end;
                    end
                end
            end
            if strcmp(fMode,'adjtrig') | strcmp(fMode,'compute') 
                load([FileBase '.' mfilename '.mat']);
                [Spec t f] = LoadSpecs(FileBase,States{w}); % t x f x ch
                Spec = exp(Spec);
                [myRes ind] = SelectPeriods(round(Res*1250/SampleRate), Period,'d',1);
                [dummy, tRes , gspk] =NearestNeighbour(t,myRes/1250,'both');
                Ph = ThPh(round(t*1250));
                [dumy fi] = min(abs(f-95));
               % chi=40;
               % c=5;
%                  atr1 = AdjustTrigAv(sq(Spec(:,fi,chi)), Ph, tRes(myClu==c), 30); 
%                  chi=7;
%                  atr2 = AdjustTrigAv(sq(Spec(:,fi,chi)), Ph, tRes(myClu==c), 30); 
                 keyboard
                for c=5%uClu'
                    if sum(myClu==c)>50
                        for fi = 1:length(f)
                            for chi =1:length(Channels)
                                myatr = AdjustTrigAv(sq(Spec(:,fi,chi)), Ph, tRes(myClu==c), 30); 
                                atr(fi,chi,:) = myatr.TrigDiffMean;
%                                  Trig = tRes(myClu==c);
%                                  bTrig = zeros(length(Ph),1);
%                                  bTrig(Trig)=1;
%                                  pc(fi,chi,c) = partialcorr(log(sq(Spec(:,fi,chi))),bTrig,log(sq(Spec(:,fi,40))),'type','Spearman');
            
                            end
                        end
                    end
                end
                 keyboard
            end

        end
        save([FileBase '.' mfilename '.mat'],'out');
        
%%%%%%%%%%%%%%%%%%%%%%%5 new method
    case 'compute_rank'
        [Res, Clu, Map ] = LoadCluRes(FileBase, Electrodes);
        [Res, Clu] = RemoveBursts(Res,Clu,10*Par.SampleRate/1000);
        Channels = RepresentChan(FileBase);
        load([FileBase '.' mfilename '.mat']);
        for w=1:length(out)
            [Spec t f] = LoadSpecs(FileBase,out(w).State); % t x f x ch
            %sSpec = exp(Spec);
            Period = load([FileBase '.sts.' out(w).State]);
            [myRes ind] = SelectPeriods(round(Res*1250/SampleRate), Period,'d',1);
             myClu = Clu(ind);
            [uClu dummy myClu] = unique(myClu);
            nClu = length(uClu);
            [dummy, tRes , gspk] =NearestNeighbour(t,myRes/1250,'both');
            %Ph = ThPh(round(t*1250));
            keyboard
           
           
        end




    case {'display','select'}
        load([FileBase '.' mfilename '.mat']);

        nSt = length(out);
        CluLoc = ClusterLocation(FileBase,Electrodes);
        [RepCh IfTetrode Info] = RepresentChan(Par);
        Map = SiliconMap(Par);
        for w=1:nSt
           figure(322121);clf
           n = size(out(w).Map,1);
           nCh = length(out(w).Channels);
           if strcmp(fMode,'select')
            out(w).Group = struct([]);
           end
           for c=1:n
               fprintf('cell %d , El = %d, Clu= %d\n',c,out(w).Map(c,1), out(w).Map(c,2));
              myind = find(CluLoc(:,1)==out(w).Map(c,1) & CluLoc(:,2)==out(w).Map(c,2));
              CellCh = CluLoc(myind,3);
              SomaChanInd = find(out(w).Channels==CellCh);
              SomaXY = Map.GridCoord(CluLoc(myind,3),:);
              UseChan = setdiff([1:nCh],SomaChanInd);
              TLag = 201;
              mat = sq(out(w).Fast.TrAvPow(:,UseChan,TLag,c)); % freq x chan
              [dummy maxi] = maxn(mat);
              MaxFi= maxi(1); mchi = maxi(2);
              MaxChi = UseChan(mchi);
              MaxCh = out(w).Channels(MaxChi);
              MaxXY = Map.GridCoord(out(w).Channels(MaxChi),:);
              while 1
                  clf
                  subplot(221);
                  mat = sq(out(w).Fast.TrAvPow(:,UseChan,TLag,c)); % freq x chan
                  % smat = conv2(gausswin(3,1),gausswin(3,1),mat,'same'); %smooth for peak detection
                  %              imagesc(out(w).Fast.f,[1:nCh-1],mat');
                  imagesc(out(w).Fast.f,[1:nCh-1],mat');
                  xlabel('Frequency'); ylabel('Channels');

                  hold on
                  Lines(out(w).Fast.f(MaxFi),[],'k');
                  Lines([],mchi,'k');

                  mappow(UseChan) = mat(MaxFi,:);
                  mappow(SomaChanInd) = NaN;
                  subplot(222);
                  MapSilicon(mappow,Par);
                  xlabel('ML'); ylabel('depth');
                  hold on
                  plot(MaxXY(1),MaxXY(2),'xk');
                  plot(SomaXY(1),SomaXY(2),'ok');

                  subplot(223)
                  %TrAvPow = freq x ch x tlag x cell
                  imagesc(out(w).Fast.tbin,out(w).Fast.f,sq(out(w).Fast.TrAvPow(:,MaxChi,:,c))); axis xy
                  xlim([-0.5 0.5]);
                  xlabel('time lag'); ylabel('Frequency');
                  hold on
                  Lines([],out(w).Fast.f(MaxFi),'k');
                  Lines(out(w).Fast.tbin(TLag),[],'k');

                  subplot(224)
                  imagesc(out(w).Fast.phbin,out(w).Fast.f,sq(out(w).Fast.Pow2Ph(:,MaxChi,c,:))); axis xy
                  xlabel('theta phase'); ylabel('Frequency');
                  Lines([],out(w).Fast.f(MaxFi),'k');

                  [x y b] = ginput(1);
                  
                  switch b
                      case 1
                          ax = WhereIClicked;
                          switch ax
                              case 1
                                  [dummy MaxFi] = min(abs(out(w).Fast.f-x));
                                  [dummy mchi] = min(abs([1:nCh-1]-y));
                              case 2
                                  xydist = dotdot(Map.GridCoord(out(w).Channels(UseChan),:),'-',[x y]);
                                  [dummy mchi] = min(sqrt(sum(xydist.^2,2)));
                              case 3
                                  [dummy MaxFi] = min(abs(out(w).Fast.f-y));
                                  [dummy TLag] = min(abs(out(w).Fast.tbin-x));
                          end
                      case 2
                          %return;
                          prteps(['sm9603_006-12_Power2Units' num2str(out(w).Map(c,:))]);
                      case 3
                          break;
                      otherwise 
                        %  keyboard;
                        if strcmp(fMode,'select')
                            out(w).Group(end+1).Freq =out(w).Fast.f(MaxFi);
                            out(w).Group(end).Tlag =out(w).Fast.tbin(TLag);
                            out(w).Group(end).Map = out(w).Map(c,:);
                            % MaxCh
                            % MaxXY
                            xyzMax = ChannelsXYZ(FileBase,MaxCh,Par);
                            xyzCell = ChannelsXYZ(FileBase,CellCh,Par);

                            out(w).Group(end).DistXY = sqrt(sum((xyzMax(1:2)-xyzCell(1:2)).^2));
                            out(w).Group(end).DistXYZ = sqrt(sum((xyzMax-xyzCell).^2));
                            out(w).Group(end).DistZ =abs(xyzMax(3)-xyzCell(3));
                            out(w).Group(end).PowMap = sq(out(w).Fast.TrAvPow(MaxFi,:,TLag,c));
                            out(w).Group(end).PowTime = sq(out(w).Fast.TrAvPow(MaxFi,MaxChi,:,c));
                            out(w).Group(end).PowSpec = sq(out(w).Fast.TrAvPow(:,MaxChi,TLag,c));
                            out(w).Group(end).xyzMax = xyzMax;
                            out(w).Group(end).xyzCell = xyzCell;
                            out(w).Group(end).CellXY = SomaXY;
                            out(w).Group(end).PowXY = MaxXY;
                            out(w).Group(end)
                        end

                  end
                  MaxChi = UseChan(mchi);
                  MaxXY = Map.GridCoord(out(w).Channels(MaxChi),:);
                  MaxCh = out(w).Channels(MaxChi);
              end %while 1 loop
              
              %waitforbuttonpress
           end %loop accross cells
           
           
        end% loop accross states

        if strcmp(fMode,'select')
            save([FileBase '.' mfilename '.mat'],'out');
        end
  
    case 'group_display'
        
        cout = CatStruct(out.Group);
        ncell = size(cout.Map,2);
        rnd = (rand(2,ncell)-0.5)*0.6;
        cxi = find(cout.PowXY(2,:)<12);
        figure(233);clf
        scatter(cout.PowXY(1,cxi)+rnd(1,cxi),cout.PowXY(2,cxi),40,cout.Freq(cxi),'filled');
        axis ij
         xlim([0.5 6.5]);ylim([0.5 16.5]);
         colorbar
        hold on
        hl = line([cout.PowXY(1,cxi)+rnd(1,cxi); cout.CellXY(1,cxi)], [cout.PowXY(2,cxi); cout.CellXY(2,cxi)]);
        set(hl,'Color','k');

        for i=1:4
            figure(233)
            in = ClusterPoints([cout.PowXY(1,:)+rnd(1,:) ; cout.PowXY(2,:)]',0);
            figure(333);
            subplot(2,2,i);
            MapSilicon(mean(unity(cout.PowMap(:,in)),2),Par);
            hold on
            hl = line([cout.PowXY(1,in)+rnd(1,in); cout.CellXY(1,in)], [cout.PowXY(2,in); cout.CellXY(2,in)]);
            set(hl,'Color','k');
        end

%         for n=1:ncell
%             
%             
%         end
%     case 'fixgroup'
%          load([FileBase '.' mfilename '.mat'],'out');
%           nSt = length(out);
%        
%         for w=1:nSt
%                 ng = length(out(w).Group);
%                 for n=1:ng
%                     c = find(out(w).Map(
%             end %loop accross cells
%         end% loop accross states
%         
% %         save([FileBase '.' mfilename '.mat'],'out');
         
end













return

