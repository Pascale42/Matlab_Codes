% function PhaseShift(FileBase, fMode, State, Channels, WinSec, Name, Ref, FreqRange)

% fMode :: 'compute' or  'display'
% State :: 'THE'... or leave empty if not needed
% Channels ::  channels for the eeg triggered representation
% WinSec :: window in seconds for the triggered average
% Name :: name it as you wish
% Ref ::  reference channel  for the trigger
%  FreqRange ::  Frequency Range to filter the signal
%
% 18 juin 2015

function PhaseShift(FileBase, fMode, State, Channels, WinSec, Name, Ref, FreqRange)


switch fMode
    
    case 'compute'
        
        Par=LoadPar([FileBase '.xml']);
        eeg=LoadBinary([FileBase '.eeg'], Channels, Par.nChannels);
        if ~isempty(State)
            STA=load([FileBase '.sts.' State]);
            eeg=SelectPeriods(eeg(:,:), STA, 'c', 1);
        end
        
        feeg = ButFilter(eeg,  2, FreqRange/(Par.lfpSampleRate/2),'bandpass');
        
        hilb = hilbert(feeg);
        ph = angle(hilb);
        
        allm = LocalMinima(eeg(:,Ref),50,0);
        
        % bound=[6.5 8.5];
        % glm = allm(log(abs(eeg(allm,Ref)))>bound(1) & log(abs(eeg(allm,Ref)))<bound(2));
        
        win=round(WinSec*Par.lfpSampleRate);
        
        % retirer les indices de bords
        a=find(allm >win,1, 'first');
        b=find(allm +win <length(feeg), 1, 'last');
        allm=allm(a:b);
        clear a b
        
        % trigered average
        nChannels=size(feeg,2);
        Av=NaN(nChannels, ((win*2)+1));
        for nch=1:nChannels
            Tot=NaN(length(glm), ((win*2)+1));
            for n=1:length(glm)
                Tot(n,:)=feeg(glm(n)-win : glm(n)+win,nch) ;
            end
            Av(nch, :)=mean(Tot,1);
        end
        
        % Faire l'axe x
        Srange = -win(1):win(end);
        Trange = linspace(-win(1),win(end),length(Srange));
        Trange=(Trange./Par.lfpSampleRate).*1000;
        clear nChannels nch Srange win WinSec Tot
        
        
        
        PhShift=NaN(size(ph,2), 1);
        for n=1:size(ph,2)
            PhShift(n)= mean(ph(allm, 1)) - mean(ph(allm, n));
        end
        
        save([FileBase '.' mfilename '.' Name '.mat'], 'Av', 'Trange', 'PhShift')
        
    case 'display'
        
        load([FileBase '.' mfilename '.' Name '.mat'], 'Av', 'Trange', 'PhShift');
        
        figure('name', ['PhaseShift - ' Name],'NumberTitle','off')
        subplot(1,4,1:3)
        imagesc(1, 1:length(PhShift), PhShift*180/pi);
        axis off; colormap('HSV'); clim([0 360]); colorbar; hold on
        PlotManyCh(Av(:,:)',Trange, 1250, 2 , 'k', 0, 0); axis off
        subplot 144
        plot(PhShift*180/pi, 1:length(PhShift),'k'); axis ij
        
end