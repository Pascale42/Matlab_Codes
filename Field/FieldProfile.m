function out = FieldProfile(FileBase,varargin)
%function FieldProfile(FileBase,Mode, FreqRange, State, SigType)
% computes power,phase profile of the LFP for various freq.
%mode can be compute/display/.. your own:)
[Mode, FreqRange, State, SigType] = DefaultArgs(varargin,{'compute',[6 9],'REM','eeg'});

eFs = 1250;
rCoef = 5;
reFs = eFs/rCoef;
MAXPER=1250*60*20;
Par = LoadPar([FileBase '.xml']);
switch Mode
    case 'compute'
        Chan = RepresentChan(FileBase);

        Period = SelectStates(FileBase,State,1250*10,0);
        dPer = diff(Period,1,2);
        %chose the period that is longest but not too much
        [maxPer maxPeri] = max(dPer);
        MyPeriod = Period(maxPeri,:);
        if maxPer>MAXPER
            MyPeriod = [MyPeriod(:,2)-MAXPER MyPeriod(:,2)];
        end
%         [sdPer si] = sort(dPer);        
%         peri =  find(sdPer<1250*60*10); %less then 10 min
%         if ~isempty(peri)sm9603_006-12
%             myper = si(peri(end));
%             MyPeriod = Period(myper,:);
%         else
%             MyPeriod = Period(1,1)+[0 eFs*60*10-1];
%         end
        if strcmp(SigType,'csd')
            if ~FileExists([FileBase '.csd'])
                 FileCSDAllCh(FileBase, Par.nChannels, Chan);               
            end
            eeg = LoadBinary([FileBase '.csd'],Chan,Par.nChannels,MyPeriod);
        else
            eeg = LoadBinary([FileBase '.eeg'],Chan,Par.nChannels,MyPeriod);
        end

        reeg = resample(double(eeg'),1,5);
        clear eeg;
        [OutArgs.y,OutArgs.f,OutArgs.ph] = mtchd(reeg,2^nextpow2(4000/rCoef),reFs,2^nextpow2(3000/rCoef),[],2);
        OutArgs.Chan = Chan;
        if nargout<1
            save([FileBase '.FieldProfile.' State '.mat'],'OutArgs');
        end
    case {'display','select'}
        load([FileBase '.FieldProfile.' State '.mat']);
        HpcChan= load([FileBase '.eegseg.par'])+1;
        HpcChan = HpcChan(1);

        Chan = OutArgs.Chan;
        HpcChani = find(Chan<=HpcChan);
        HpcChani  = HpcChani(end);
        CxChani = find(Chan<=HpcChan);
        nCh = length(Chan);
        for i=1:nCh
            pow(:,i) = (sq(OutArgs.y(:,i,i)));
        end
        myf = find(OutArgs.f>FreqRange(1) & OutArgs.f<FreqRange(2));
       
     %   [dummy maxch] = max(mypow);
        maxch =HpcChani;
        ph = sq(OutArgs.ph(:,maxch,:));
        
        coh = sq(OutArgs.y(:,maxch,:));
        coh(:,maxch) = ones(size(coh,1),1);
        mycoh = mean(coh(myf,:));
        myph = circmean(ph(myf,:))*180/pi;
        
        mypow = mean(pow(myf,:));
        %mypow = mypow/mypow(maxch);
        
        nmyf = setdiff(find(OutArgs.f>1 & OutArgs.f<14),myf);
        powsnr = mean(pow(myf,:))./mean(pow(nmyf,:));
        
        mycoh(HpcChani) = NaN;
        myph(HpcChani) = NaN;
        Map = SiliconMap(Par);
        
        switch Mode
            case {'display','select'}
                figure(9983);clf
                subplot(231)
                imagesc(OutArgs.f,[1:nCh],pow');
                title('power');
                xlabel('freq. (Hz)');
                set(gca,'YTick',[1:nCh]);
                set(gca,'YTickLabel',num2str(Chan(:)));
                %xlim([FreqRange]);
                colorbar

                subplot(232)
                imagesc(OutArgs.f,[1:nCh],coh');
                title(['coherence rel to ch' num2str(HpcChan) ', rad']);
                xlabel('freq. (Hz)');
                set(gca,'YTickLabel','');
                set(gca,'YTick',[1:nCh]);
                set(gca,'YTickLabel',num2str(Chan(:)));

                %xlim([FreqRange]);
                colorbar

                subplot(233)
                imagesc(OutArgs.f,[1:nCh],ph'*180/pi);
                title(['phase shift rel to ch  ' num2str(HpcChan) ', rad']);
                xlabel('freq. (Hz)');
                set(gca,'YTickLabel','');
                %xlim([FreqRange]);
                set(gca,'YTick',[1:nCh]);
                set(gca,'YTickLabel',num2str(Chan(:)));
                caxis([-180 180]);
                colorbar

                subplot(234)
                MapSilicon(mypow,OutArgs.Chan,Par,[],1,[],[],[],1);
                hold on
                plot(Map.GridCoord(HpcChan,1),Map.GridCoord(HpcChan,2),'ko','MarkerSize',10);

                subplot(235)
                MapSilicon(mycoh,OutArgs.Chan,Par,[],1,[],[],[],1);
                hold on
                plot(Map.GridCoord(HpcChan,1),Map.GridCoord(HpcChan,2),'ko','MarkerSize',10);
               
                subplot(236)
                MapSilicon(myph,OutArgs.Chan,Par,[],1,[-180 180],1,[],1);
                hold on
                plot(Map.GridCoord(HpcChan,1),Map.GridCoord(HpcChan,2),'ko','MarkerSize',10);
                caxis([-60 60]);
                
                if strcmp(Mode,'select')
                    ChMap = Map.GridCoord(OutArgs.Chan,:);
                    [xch ych b] = ginput(2);
                    for l=1:2
                        [dd chi(l)] = min(abs(ChMap(:,1)-xch(l)).^2+abs(ChMap(:,2)-ych(l)).^2);
                    end

                    chi = sort(chi);
                    OutArgs.GroupChanInd = chi;
                    OutArgs.Hpc.Pow = mypow(HpcChani);
                    OutArgs.Hpc.PowSnr = powsnr(HpcChani);

                    OutArgs.Cx.Pow = mypow(chi(1):chi(2));
                    OutArgs.Cx.PowSnr = powsnr(chi(1):chi(2));
                    OutArgs.Cx.Coh = mycoh(chi(1):chi(2));
                    OutArgs.Cx.Ph =  myph(chi(1):chi(2));
                    save([FileBase '.FieldProfile.' State '.mat'],'OutArgs');                    
                end

        end
        
    case 'group'
          load([FileBase '.FieldProfile.' State '.mat']);
       
          x =[1:length(OutArgs.Cx.Pow)];
          bcoh = robustfit(x,OutArgs.Cx.Coh);
          pcoh = polyfit(x,OutArgs.Cx.Coh,2);
          out.CohSlope = bcoh(2);
          out.CohPoly = pcoh;
          bpow = robustfit(x,sqrt(OutArgs.Cx.Pow));
          ppow = polyfit(x,sqrt(OutArgs.Cx.Pow),2);
          out.PowSlope = bpow(2);
          out.PowPoly = ppow;
          out.MeanPh = 180*circmean(OutArgs.Cx.Ph*pi/180)/pi;
          out.HpcPow = OutArgs.Hpc.Pow;
          ChDist = diff(OutArgs.Chan(OutArgs.GroupChanInd));
          out.PowDrop = (sqrt(OutArgs.Cx.Pow(end))-sqrt(OutArgs.Cx.Pow(1)))/ChDist/max(sqrt(OutArgs.Cx.Pow));
          out.CohDrop = (sqrt(OutArgs.Cx.Coh(end))-sqrt(OutArgs.Cx.Coh(1)))/ChDist;
          if 1
              figure(3232);clf
              subplot(221)
              plot(sqrt(OutArgs.Cx.Pow),'o-');
              hold on
              plot(x,bpow(1)+bpow(2)*x,'g');
              plot(x,ppow(3)+ppow(2)*x+ppow(1)*x.^2,'m');
              Lines([],sqrt(OutArgs.Hpc.Pow),'r','--');
              title('theta power');

              subplot(222)
              plot(OutArgs.Cx.PowSnr);
              hold on
              Lines([],OutArgs.Hpc.PowSnr,'r','--');
              title('theta power snr');

              subplot(223)
              plot(OutArgs.Cx.Coh,'o-');
              hold on
              plot(x,bcoh(1)+bcoh(2)*x,'g');
              plot(x,pcoh(3)+pcoh(2)*x+pcoh(1)*x.^2,'m');
              title('coherence');

              subplot(224)
              plot(OutArgs.Cx.Ph);
              title('phase');
          end
    
    case 'group2'
          load([FileBase '.FieldProfile.' State '.mat']);
          HpcChan= load([FileBase '.eegseg.par'])+1;
          HpcChan = HpcChan(1);

         Chan = OutArgs.Chan;
         [dummy HpcChani] = min(abs(Chan-HpcChan));
        
         [dummy fi] =max( OutArgs.y(:,HpcChani,HpcChani));
         out.f = OutArgs.f(fi);
         out.PowCa1 = sqrt((OutArgs.y(fi,HpcChani,HpcChani)));
         out.PowCa3 = sqrt((OutArgs.y(fi,end,end)));
         out.CohCa1Ca3 = ((OutArgs.y(fi,HpcChani,end)));
         out.PhCa1Ca3 = ((OutArgs.ph(fi,HpcChani,end)));
end

