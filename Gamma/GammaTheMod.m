% function GammaTheMod(eeg,RefCh,nChannels,FreqR,W, Mode, fMode)

function GammaTheMod(eeg,RefCh,nChannels,FreqR,W,Mode,fMode)

switch Mode

  %% Modulation of gamma POWER
  case 'power'

    switch fMode
      case 'compute'

        SR=1250;
        WinSec=0.05; % let's narrow the window to have more precise points for phase
        WinSp = 2^round(log2(WinSec*SR)); % 0.082 s ; so 0.041s overlap
        nFFT = 2*WinSp; % 1.6 s
        if W==1
          eeg=WhitenSignal(eeg,[],1);
        else
        end
%         [y,f,t]= mtchglong_K(eeg, nFFT,SR,WinSp,[],3,'linear',[], FreqR);
        [y,f,t]= mtchglong(eeg, nFFT,SR,WinSp,[],3,'linear',[], FreqR);

        pow=zeros(length(t),length(f),nChannels);
        for n=1:nChannels
          pow(:,:,n)=20*log10(abs(y(:,:,n,n))+eps);
        end

        % get mean power over all the frequency band 'FreqR'
        powm=squeeze(mean(pow,2));


        % get instantaneous phase regarding RefCh for each power-time binls
        feeg= ButFilter(eeg(:,RefCh), 2, [2 5]/(1250/2),'bandpass');
        hilb= hilbert(feeg);
        ph = angle(hilb);
        res=round(t.*1250);
        phres=ph(res);

        % Sort and bin
        [phsort id]=sort(phres,1);
        phpow=powm(id,:);
        phpow=[phsort phpow];
        phpow(:,1)=phpow(:,1).*180/pi;


        % mean power over 5 degres bins
        int=[[-180:5:175]' [-175:5:180]'];
        powbin=zeros(72,nChannels);
        for n=1:size(int,1)
          bin=find(phpow(:,1)>= int(n,1) & phpow(:,1) <int(n,2));
          for ch=1:nChannels
            powbin(n,ch)=mean(phpow(bin,ch+1));
          end
        end

        if FreqR(2)>80; oscill='vhf'; else; oscill='gam'; end
        save([mfilename Mode '_chr' num2str(RefCh) oscill '.mat'],'t','f','y','feeg','phpow','powbin','nChannels')


      case 'display'

        if FreqR(2)>80; oscill='vhf'; else; oscill='gam'; end
        load([mfilename Mode '_chr' num2str(RefCh) oscill '.mat'])


%         figure('Name', 'Theta modul of gamma power in time','NumberTitle','off')
%         subplot(211)
%         imagesc(t,f,20*log10(abs(y(:,:,RefCh,RefCh))+eps)');
%         axis xy; colormap(jet); title(['Power specgram ' num2str(RefCh)]); xlabel('Time'); ylabel('Frequency')
%         subplot(212)
%         plot(feeg); axis tight
% 
%         figure('Name','Distrib power as fct of phase for each channel','NumberTitle','off');
%         for i=1:nChannels; subplot(6,6,i); plot(phpow(:,1),phpow(:,i+1),'o');end
%         copyax('on')

        % figure('Name','theta modulation of gamma power','NumberTitle','off');
        % imagesc([phpow(:,1);(phpow(:,1)+360) ],[1:nChannels],[phpow(:,2:nChannels+1);phpow(:,2:nChannels+1)]')
        int=[[-180:5:175]' [-175:5:180]'];

        figure('Name','theta modulation of gamma power','NumberTitle','off'); imagesc([int(:,1);[182.5:5:537.5]'],[1:nChannels],[powbin;powbin]')
        set(gca,'XTick',[-180:90:540]); colorbar
        ylabel('Channels'); xlabel('Theta Phase') % super!

    end
    
    
    
      %% Modulation of gamma POWER
  
    
    
      %% Modulation of gamma POWER freq for one channel
    case 'powerdc'

    switch fMode
      case 'compute'

        SR=1250;
        WinSec=0.05; % let's narrow the window to have more precise points for phase
        WinSp = 2^round(log2(WinSec*SR)); % 0.082 s ; so 0.041s overlap
        nFFT = 2*WinSp; % 1.6 s
        if W==1
          eeg=WhitenSignal(eeg,[],1);
        else
        end
        
        
        
         % get instantaneous phase regarding RefCh for each power-time binls
        feeg= ButFilter(eeg(:,RefCh), 2, [2 5]/(1250/2),'bandpass');
        hilb= hilbert(feeg);
        ph = angle(hilb);
        res=round(t.*1250);
        phres=ph(res);
        
        
%         [y,f,t]= mtchglong_K(eeg, nFFT,SR,WinSp,[],3,'linear',[], FreqR);

        Eeg=eegt(1:1250*1200,[3 6 9 14]);%###########
        [Y,F,T]= mtchglong(Eeg, nFFT,SR,WinSp,[],3,'linear',[], [20 100]);

        nChannels=4;%###########
        
        pow=zeros(length(T),length(F),nChannels);
        for n=1:nChannels
          pow(:,:,n)=20*log10(abs(Y(:,:,n,n))+eps);
        end

        % get all the power for one channel
        
        pow1=pow(:,:,1);
        pow2=pow(:,:,2);
        pow3=pow(:,:,3);
        pow4=pow(:,:,4);

       

        % Sort and bin
        [phsort id]=sort(phres,1);
        
        phpow1=pow1(id,:);
        phpow1=[phsort phpow1];
        phpow1(:,1)=phpow1(:,1).*180/pi;
        
        phpow2=pow2(id,:);
        phpow2=[phsort phpow2];
        phpow2(:,1)=phpow2(:,1).*180/pi;
        
        phpow3=pow3(id,:);
        phpow3=[phsort phpow3];
        phpow3(:,1)=phpow3(:,1).*180/pi;

        phpow4=pow4(id,:);
        phpow4=[phsort phpow4];
        phpow4(:,1)=phpow4(:,1).*180/pi;
        
        nChannels = size(pow1,2);

        
        
        % mean power over 5 degres bins
        int=[[-180:5:175]' [-175:5:180]'];
        powbin1=zeros(72,nChannels);
        powbin2=zeros(72,nChannels);
        powbin3=zeros(72,nChannels);
        powbin4=zeros(72,nChannels);        
        
        for n=1:size(int,1)
          bin1=find(phpow1(:,1)>= int(n,1) & phpow1(:,1) <int(n,2));
          bin2=find(phpow2(:,1)>= int(n,1) & phpow2(:,1) <int(n,2));
          bin3=find(phpow3(:,1)>= int(n,1) & phpow3(:,1) <int(n,2));
          bin4=find(phpow4(:,1)>= int(n,1) & phpow4(:,1) <int(n,2));
          for ch=1:nChannels
            powbin1(n,ch)=mean(phpow1(bin1,ch+1));
            powbin2(n,ch)=mean(phpow2(bin2,ch+1));
            powbin3(n,ch)=mean(phpow3(bin3,ch+1));
            powbin4(n,ch)=mean(phpow4(bin4,ch+1));
          end
        end

      save([FileBase '.' mfilename 'GamPowTheMod_FreqLayers-20-100Hz.mat'],'T','F','Y','phpow1','powbin1','phpow2','powbin2','phpow3','powbin3','phpow4','powbin4')


      case 'display'

     
        load([FileBase '.' mfilename 'GamPowTheMod_FreqLayers.mat']);

        int=[[-180:5:175]' [-175:5:180]'];

        figure;
        subplot(4,1,1)
        imagesc([int(:,1);[182.5:5:537.5]'],F,[powbin1;powbin1]');
        set(gca,'XTick',[-180:90:540]); colorbar; axis xy;title('Layers1 - 2 up - 2 bottom - 3'); clim([60 90])
        subplot(4,1,2);
        imagesc([int(:,1);[182.5:5:537.5]'],F,[powbin2;powbin2]');
        set(gca,'XTick',[-180:90:540]); colorbar; axis xy; clim([60 90])
        subplot(4,1,3);
        imagesc([int(:,1);[182.5:5:537.5]'],F,[powbin3;powbin3]');
        set(gca,'XTick',[-180:90:540]); colorbar; axis xy; clim([60 90])
        subplot(4,1,4);
        imagesc([int(:,1);[182.5:5:537.5]'],F,[powbin4;powbin4]');
        set(gca,'XTick',[-180:90:540]); colorbar; axis xy; clim([60 90])
        
                ylabel('Channels'); xlabel('Theta Phase') % super!
    end
    
    
    
    
%% Modulation of gamma POWER NORMALIZED
  
    case 'pownorm'

    switch fMode
      case 'compute'

         if FreqR(2)>80; oscill='vhf'; else; oscill='gam'; end
        
        if FileExists([mfilename 'power_chr' num2str(RefCh) oscill '.mat'])
          load([mfilename 'power_chr' num2str(RefCh) oscill '.mat'],'t','f','y','feeg')
        else
        
        
        SR=1250;
        WinSec=0.05; % let's narrow the window to have more precise points for phase
        WinSp = 2^round(log2(WinSec*SR)); % 0.082 s ; so 0.041s overlap
        nFFT = 2*WinSp; % 1.6 s
        FreqR=[20 55];

        [y,f,t]= mtchglong(eeg, nFFT,SR,WinSp,[],3,'linear',[], FreqR);
         feeg= ButFilter(eeg(:,RefCh), 2, [2 5]/(1250/2),'bandpass');
        end

        pow=zeros(length(t),length(f),nChannels);
        for n=1:nChannels
          pow(:,:,n)=20*log10(abs(y(:,:,n,n))+eps);
        end

        % get mean power over all the frequency band 'FreqR'
        powm=squeeze(mean(pow,2));

        
        % normalize power for layer 1 gamma power
        powmax=zeros(size(powm,1),size(powm,2));
        for n=1:nChannels
          powmax(:,n)=powm(:,n)./ max(powm(:,n));
        end

        % get instantaneous phase regarding RefCh for each power-time bin
       
        hilb= hilbert(feeg);
        ph = angle(hilb);
        res=round(t.*1250);
        phres=ph(res);

        % Sort and bin
        [phsort id]=sort(phres,1);
        phpowmax=powmax(id,:);
        phpowmax=[phsort phpowmax];
        phpowmax(:,1)=phpowmax(:,1).*180/pi;


        % mean power over 5 degres bins
        int=[[-180:5:175]' [-175:5:180]'];
        powmaxbin=zeros(72,nChannels);
        for n=1:size(int,1)
          bin=find(phpowmax(:,1)>= int(n,1) & phpowmax(:,1) <int(n,2));
          for ch=1:nChannels
            powmaxbin(n,ch)=mean(phpowmax(bin,ch+1));
          end
        end

        if FreqR(2)>80; oscill='vhf'; else; oscill='gam'; end
        save([mfilename Mode '_chr' num2str(RefCh) oscill '.mat'],'t','f','y','phpowmax','powmaxbin','nChannels')


      case 'display'

        if FreqR(2)>80; oscill='vhf'; else; oscill='gam'; end
        load([mfilename Mode '_chr' num2str(RefCh) oscill '.mat'])


%         figure('Name', 'Theta modul of gamma power in time','NumberTitle','off')
%         subplot(211)

%         imagesc(t,f,20*log10(abs(y(:,:,RefCh,RefCh))+eps)');
%         axis xy; colormap(jet); title(['Power specgram ' num2str(RefCh)]); xlabel('Time'); ylabel('Frequency')
%         subplot(212)
%         plot(feeg); axis tight
% 
%         figure('Name','Distrib power as fct of phase for each channel','NumberTitle','off');
%         for i=1:nChannels; subplot(6,6,i); plot(phpow(:,1),phpow(:,i+1),'o');end
%         copyax('on')

        % figure('Name','theta modulation of gamma power','NumberTitle','off');
        % imagesc([phpow(:,1);(phpow(:,1)+360) ],[1:nChannels],[phpow(:,2:nChannels+1);phpow(:,2:nChannels+1)]')
        int=[[-180:5:175]' [-175:5:180]'];

        figure('Name','theta modulation of gamma power','NumberTitle','off'); imagesc([int(:,1);[182.5:5:537.5]'],[1:nChannels],[powmaxbin;powmaxbin]')
        set(gca,'XTick',[-180:90:540]); colorbar
        ylabel('Channels'); xlabel('Theta Phase') % super!

    end

%% Modulation of gamma COHERENCE

  case 'coherence'

    switch fMode
      case 'compute'

        if FreqR(2)>80; oscill='vhf'; else; oscill='gam'; end
        
        if FileExists([mfilename 'power_chr' num2str(RefCh) oscill '.mat'])
          load([mfilename Mode 'power_chr' num2str(RefCh) oscill '.mat'],'t','f','y','feeg')
        else

          SR=1250;
          WinSec=0.05; % let's narrow the window to have more precise points for phase
          WinSp = 2^round(log2(WinSec*SR)); % 0.082 s ; so 0.041s overlap
          nFFT = 2*WinSp; % 1.6 s
          FreqR=[20 55];

          [y,f,t]= mtchglong_K(eeg, nFFT,SR,WinSp,[],3,'linear',[], FreqR);
          feeg= ButFilter(eeg(:,RefCh), 2, [2 5]/(1250/2),'bandpass');
        end
        
        
        coh=zeros(length(t),length(f),nChannels);
        for n=1:nChannels
          coh(:,:,n)=y(:,:,n,RefCh);
        end

        % get mean coherence over all the frequency band 'FreqR'
        cohm=squeeze(mean(coh,2));
        cohm(:,RefCh)=1;

        % get instantaneous phase regarding RefCh for each power-time bin
        
        hilb= hilbert(feeg);
        ph = angle(hilb);
        res=round(t.*1250);
        phres=ph(res);
      
        % Sort and bin
        [phsort id]=sort(phres,1);
        phcoh=cohm(id,:);
        phcoh=[phsort phcoh];
        phcoh(:,1)=phcoh(:,1).*180/pi;


        % mean power over 5 degres bins
        int=[[-180:5:175]' [-175:5:180]'];
        cohbin=zeros(72,nChannels);
        for n=1:size(int,1)
          bin=find(phcoh(:,1)>= int(n,1) & phcoh(:,1) <int(n,2));
          for ch=1:nChannels
            cohbin(n,ch)=mean(phcoh(bin,ch+1));
          end
        end



        save([mfilename Mode '_chr' num2str(RefCh) oscill '.mat'],'t','f','y','feeg','phcoh','cohbin','nChannels')


      case 'display'

        if FreqR(2)>80; oscill='vhf'; else; oscill='gam'; end
        load([mfilename Mode '_chr' num2str(RefCh) oscill '.mat'],'t','f','y','phcoh','cohbin')


% 
%         figure('Name','Distrib coherence as fct of phase for each channel','NumberTitle','off');
%         for i=1:nChannels; subplot(6,6,i); plot(phcoh(:,1),phcoh(:,i+1),'o');end
%         copyax('on')

        % figure('Name','theta modulation of gamma power','NumberTitle','off');
        % imagesc([phpow(:,1);(phpow(:,1)+360) ],[1:nChannels],[phpow(:,2:nChannels+1);phpow(:,2:nChannels+1)]')
        int=[[-180:5:175]' [-175:5:180]'];

        figure('Name','theta modulation of gamma coherence','NumberTitle','off'); imagesc([int(:,1);[182.5:5:537.5]'],[1:nChannels],[cohbin;cohbin]')
        set(gca,'XTick',[-180:90:540]); colorbar
        ylabel('Channels'); xlabel('Theta Phase') % super!

    end
    
%% Modulation of gamma COHERENCE NORMALIZED

  case 'cohnorm'

    switch fMode
      case 'compute'

        if FreqR(2)>80; oscill='vhf'; else; oscill='gam'; end
        
        if FileExists([mfilename Mode '_chr' num2str(RefCh) oscill '.mat'])
          load([mfilename Mode '_chr' num2str(RefCh) oscill '.mat'],'t','f','y','feeg')
        else

          SR=1250;
          WinSec=0.05; % let's narrow the window to have more precise points for phase
          WinSp = 2^round(log2(WinSec*SR)); % 0.082 s ; so 0.041s overlap
          nFFT = 2*WinSp; % 1.6 s
          FreqR=[20 55];

          [y,f,t]= mtchglong_K(eeg, nFFT,SR,WinSp,[],3,'linear',[], FreqR);
          feeg= ButFilter(eeg(:,RefCh), 2, [2 5]/(1250/2),'bandpass');
        end
        
        
        coh=zeros(length(t),length(f),nChannels);
        for n=1:nChannels
          coh(:,:,n)=y(:,:,n,RefCh);
        end

        % get mean coherence over all the frequency band 'FreqR'
        cohm=squeeze(mean(coh,2));
        cohm(:,RefCh)=1;

        % get instantaneous phase regarding RefCh for each power-time bin
        
        hilb= hilbert(feeg);
        ph = angle(hilb);
        res=round(t.*1250);
        phres=ph(res);
      
        % Sort and bin
        [phsort id]=sort(phres,1);
        phcoh=cohm(id,:);
        phcoh=[phsort phcoh];
        phcoh(:,1)=phcoh(:,1).*180/pi;


        % mean power over 5 degres bins
        int=[[-180:5:175]' [-175:5:180]'];
        cohbin=zeros(72,nChannels);
        for n=1:size(int,1)
          bin=find(phcoh(:,1)>= int(n,1) & phcoh(:,1) <int(n,2));
          for ch=1:nChannels
            cohbin(n,ch)=mean(phcoh(bin,ch+1));
          end
        end



        save([mfilename Mode '_chr' num2str(RefCh) oscill '.mat'],'t','f','y','feeg','phcoh','cohbin','nChannels')


      case 'display'

        if FreqR(2)>80; oscill='vhf'; else; oscill='gam'; end
        load([mfilename Mode '_chr' num2str(RefCh) oscill '.mat'],'t','f','y','phcoh','cohbin')


% 
%         figure('Name','Distrib coherence as fct of phase for each channel','NumberTitle','off');
%         for i=1:nChannels; subplot(6,6,i); plot(phcoh(:,1),phcoh(:,i+1),'o');end
%         copyax('on')

        % figure('Name','theta modulation of gamma power','NumberTitle','off');
        % imagesc([phpow(:,1);(phpow(:,1)+360) ],[1:nChannels],[phpow(:,2:nChannels+1);phpow(:,2:nChannels+1)]')
        int=[[-180:5:175]' [-175:5:180]'];

        figure('Name','theta modulation of gamma coherence','NumberTitle','off'); imagesc([int(:,1);[182.5:5:537.5]'],[1:nChannels],[cohbin;cohbin]')
        set(gca,'XTick',[-180:90:540]); colorbar
        ylabel('Channels'); xlabel('Theta Phase') % super!

    end
    
%% Modulation of gamma BURSTS

  case 'counts'

    switch fMode
      case 'compute'
        thefeeg = ButFilter(eeg(:,RefCh), 2, [2 5]/(1250/2),'bandpass');
        gamfeeg= ButFilter(eeg, 2, FreqR/(1250/2),'bandpass');

        hilbr = hilbert(thefeeg);
        phr = angle(hilbr);

        phhist=zeros(nChannels,72); % 72 bins for phase


        for ch=1:nChannels
          sigma=std(gamfeeg(:,ch));
          if ch>=1 & ch<=3; coef=0.5; else; coef=2;end
          x=find(gamfeeg(:,ch)>= coef*sigma & gamfeeg(:,ch)>= -coef*sigma);
          gamtheph{ch}=phr(x).*180/pi;
          [phhist(ch,:) phbins]=hist(gamtheph{ch},[-177.5:5:177.5]);
        end
        phhist=phhist';

        if FreqR(2)>80; oscill='vhf'; else; oscill='gam'; end

        save([mfilename Mode '_chr' num2str(RefCh) oscill '.mat'],'phhist', 'phbins','gamtheph')

        %
        % % just for test gamma in time
        % teeg=eeg(1:12500,:);
        % SR=1250;
        % WinSec=0.05; % let's narrow the window to have more precise points for phase
        % WinSp = 2^round(log2(WinSec*SR)); % 0.82 s ; so 0.41s overlap
        % nFFT = 2*WinSp; % 1.6 s
        % FR=[20 55];
        % [y,f,t]= mtchglong_K(teeg, nFFT,SR,WinSp,[],3,'linear',[], FR);
        %
        % figure
        %   imagesc(t,f,20*log10(abs(y(:,:,1,1))+eps)');
        %   axis xy; colormap(jet); title(['Power specgram ' num2str(n)]); xlabel('Time'); ylabel('Frequency')
        %
        % meangam=zeros(size(y,1),31);
        % for n=1:31
        %   meangam(:,n)=mean(y(:,:,n,n),2);
        % end
        %
        % figure
        %   subplot(211)
        %   imagesc(t,[1:31],20*log10(abs(meangam)+eps)');
        % colormap(jet);xlabel('Time'); ylabel('Frequency')
        %   subplot(212)
        %   plot(feeg(1:12500,1)); axis tight
        %
        %
        %   % just for test vhf in time
        % teeg=eeg(1:12500,:);
        % SR=1250;
        % WinSec=0.05; % let's narrow the window to have more precise points for phase
        % WinSp = 2^round(log2(WinSec*SR)); % 0.82 s ; so 0.41s overlap
        % nFFT = 2*WinSp; % 1.6 s
        % FR=[250 400];
        % [yv,fv,tv]= mtchglong_K(teeg, nFFT,SR,WinSp,[],3,'linear',[], FR);
        %
        % figure
        %   imagesc(tv,fv,20*log10(abs(yv(:,:,1,1))+eps)');
        %   axis xy; colormap(jet); title(['Power specgram ' num2str(14)]); xlabel('Time'); ylabel('Frequency')
        %
        % meanvhf=zeros(389,31);
        % for n=1:31
        %   meanvhf(:,n)=mean(yv(:,:,n,n),2);
        % end
        %
        % figure
        %   subplot(211)
        %   imagesc(tv,[1:31],20*log10(abs(meanvhf)+eps)');
        % colormap(jet);xlabel('Time'); ylabel('Frequency')
        %   subplot(212)
        %   plot(feeg(1:12500,1)); axis tight

      case 'display'

        if FreqR(2)>80; oscill='vhf'; else; oscill='gam'; end
        load([mfilename Mode '_chr' num2str(RefCh) oscill '.mat']);

        if FreqR(2)>80;  figure('Name', 'Theta modulation of high gamma power');
        else; figure('Name', 'Theta modulation of Gamma power'); end
        imagesc([phbins [182.5:5:537.5]],[1:nChannels],unity([phhist;phhist])');
        set(gca,'XTick',[-180:90:540])
        ylabel('Channels'); xlabel('Theta Phase')

    end
end

%
% figure;for ch=1:31; subplot(121); hist(gamtheph{ch},60); subplot(122);rose(gamtheph{ch},60); waitforbuttonpress; end
% figure;for ch=1:31; hist(gamtheph{ch},60);  waitforbuttonpress; end
% figure;for ch=1:31;  subplot(6,6,ch); hist(gamtheph{ch},60);  axis tight; end
% figure;for ch=1:31;  subplot(6,6,ch); imagesc(phbins,1,phhist(ch,:)');  axis tight; end
%