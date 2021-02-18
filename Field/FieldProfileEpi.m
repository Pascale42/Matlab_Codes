function out = FieldProfileEpi(FileBase,varargin)
%function FieldProfileEpi(FileBase,Mode, FreqRange, State, SigType)
% computes power,phase profile of the LFP for various freq.
%mode can be compute/display/.. your own:)
[Mode, FreqRange, State, SigType] = DefaultArgs(varargin,{'compute',[6 9],'REM','eeg'});

Par = LoadPar([FileBase '.xml']);
eFs = Par.lfpSampleRate;
if eFs==1250
   rCoef = 5;
else
    rCoef =1;
end
reFs = eFs/rCoef;
MAXPER=Par.lfpSampleRate*60*20;

switch Mode
    case 'compute'
        Chan = RepresentChan(FileBase);

        Period = SelectStates(FileBase,State,eFs*10,0);
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

        reeg = resample(double(eeg'),1,rCoef);
        clear eeg;
        [OutArgs.y,OutArgs.f,OutArgs.ph] = mtchd(reeg,2^nextpow2(4000*eFs/1250/rCoef),reFs,2^nextpow2(3000*eFs/1250/rCoef),[],2);
        OutArgs.Chan = Chan;
        if nargout<1
            save([FileBase '.FieldProfileEpi.' State '.mat'],'OutArgs');
        end
    case {'display','group','falloff'}
        load([FileBase '.' mfilename '.' State '.mat']);
        HpcChan= load([FileBase '.eegseg.par'])+1;
        HpcChan = HpcChan(1);

        Chan = OutArgs.Chan;
        HpcChani = find(Chan<=HpcChan);
        HpcChani  = HpcChani(end);
        CxChani = find(Chan<=HpcChan);
        nCh = length(Chan);
        for i=1:nCh
            pow(:,i) = sqrt(sq(OutArgs.y(:,i,i)));
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
%        mypow = mypow/mypow(maxch);
        
        %nmyf = setdiff([1:length(OutArgs.f)],myf);
%        nmyf = 
 %       powsnr = mean(pow(myf,:))./mean(pow(nmyf,:));
        
        mycoh(HpcChani) = NaN;
        myph(HpcChani) = NaN;
        Map = SiliconMap(Par);
        
        switch Mode
            case {'display','group'}
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
                imagesc(OutArgs.f,[1:nCh],ph'*180/pi');
                title(['phase shift rel to ch  ' num2str(HpcChan) ', rad']);
                xlabel('freq. (Hz)');
                set(gca,'YTickLabel','');
                %xlim([FreqRange]);
                set(gca,'YTick',[1:nCh]);
                set(gca,'YTickLabel',num2str(Chan(:)));
                caxis([-180 180]);
                colorbar

                EpiChani = find(OutArgs.Chan<=32);
                load([FileBase '.EpiMap.mat']);
          %      EpiMap = EpiduralMap('l',90);
                subplot(234)
                ImageEpiduralMap(mypow(EpiChani),EpiMap,'imagesc',OutArgs.Chan(EpiChani));
%                MapSilicon(mypow,OutArgs.Chan,Par,[],1);
                hold on
                plot(Map.GridCoord(HpcChan,1),Map.GridCoord(HpcChan,2),'ko','MarkerSize',10);

                subplot(235)
%                MapSilicon(mycoh,OutArgs.Chan,Par,[],1);
                ImageEpiduralMap(mycoh(EpiChani),EpiMap,'imagesc',OutArgs.Chan(EpiChani));
                hold on
                plot(Map.GridCoord(HpcChan,1),Map.GridCoord(HpcChan,2),'ko','MarkerSize',10);
               
                subplot(236)
               ImageEpiduralMap(myph(EpiChani),EpiMap,'imagesc',OutArgs.Chan(EpiChani));
%                MapSilicon(myph,OutArgs.Chan,Par,[],1,[-180 180],1);
                hold on
                plot(Map.GridCoord(HpcChan,1),Map.GridCoord(HpcChan,2),'ko','MarkerSize',10);
                caxis([-60 60]);
                
                if strcmp(Mode,'group')
%                     ChMap = Map.GridCoord(OutArgs.Chan,:);
%                     [xch ych b] = ginput(2);
%                     for l=1:2
%                         [dd chi(l)] = min(abs(ChMap(:,1)-xch(l)).^2+abs(ChMap(:,2)-ych(l)).^2);
%                     end
% 
%                     chi = sort(chi);
%                     OutArgs.GroupChanInd = chi;
%                     OutArgs.Hpc.Pow = mypow(HpcChani);
%                     OutArgs.Hpc.PowSnr = powsnr(HpcChani);
% 
%                     OutArgs.Cx.Pow = mypow(chi(1):chi(2));
%                     OutArgs.Cx.PowSnr = powsnr(chi(1):chi(2));
%                     OutArgs.Cx.Coh = mycoh(chi(1):chi(2));
%                     OutArgs.Cx.Ph =  myph(chi(1):chi(2));
                    
                end
                save([FileBase '.FieldProfile.' State '.mat'],'OutArgs');
            case 'falloff'
                    EpiChani = find(OutArgs.Chan<=32);
                    load([FileBase '.EpiMap.mat']);
                    Coords = EpiMap.Coord(OutArgs.Chan(EpiChani),:);
                    [maxamp maxch] = max(mypow(EpiChani));
                    Dist = sqrt(sum(dotdot(Coords(EpiChani,:),'-',Coords(maxch,:)).^2,2));

                    figure(3232311);clf
                    subplot(311)
                    PowDrop = 100*(1-mypow(EpiChani)'/mypow(HpcChani));
                    plot(Dist',PowDrop,'.');
                    
                    subplot(312)
                    CohDrop = 100*(1-mycoh(EpiChani));
                    plot(Dist',CohDrop','.')
                    
                    out.MeanCohSurf = mean(CohDrop(Dist<2));
                    out.MeanPowSurf = mean(PowDrop(Dist<2));
                    
                    b = robustfit(Dist,PowDrop);
                    out.SlopePowDrop = b(2);
                    
                    b = robustfit(Dist,CohDrop);
                    out.SlopeCohDrop = b(2);
                    
                    subplot(313)
                    PhDrop = myph(EpiChani);
                    plot(Dist',PhDrop','.')
                    keyboard
                    
        end
        
    case 'group'
          load([FileBase '.FieldProfile.' State '.mat']);
          
          if 0
          figure(3232);clf
          subplot(221)
          plot(OutArgs.Cx.Pow);
          hold on
          Lines([],OutArgs.Hpc.Pow,'r','--');
          title('theta power');
          
          subplot(222)
          plot(OutArgs.Cx.PowSnr);
          hold on
          Lines([],OutArgs.Hpc.PowSnr,'r','--');
          title('theta power snr');
          
          subplot(223)
          plot(OutArgs.Cx.Coh);
          title('coherence');
        
          subplot(224)
          plot(OutArgs.Cx.Ph);
          title('phase');
          end
          x =[1:length(OutArgs.Cx.Pow)];
          b = robustfit(x,OutArgs.Cx.Coh);
          out.CohSlope = b(2);
          b = robustfit(x,OutArgs.Cx.Pow);
          out.PowSlope = b(2);
          out.MeanPh = 180*circmean(OutArgs.Cx.Ph*pi/180)/pi;
          
end

