% function  TriggerCSD(FileBase, fMode, State, Channels, WinSec, Name, Ref, FreqRange)

% fMode :: 'compOscill', 'compExtTrig',  'display'
%           'compOscill' fMode :: uses an oscillatory activity to get its triggered CSD
%           'compExtTrig' fMode :: uses an external trigger to compute the CSD of the eeg
% State :: 'THE'... or leave empty if not needed
% Channels ::  channels for CSD
% WinSec :: window in seconds for the triggered average
% Name :: name it as you wish
% Ref ::  reference channel  for the trigger ('compOscill' fMode) OR
%            your external triggers (times in samples, @1250Hz) ('compExtTrig' fMode)
%  FreqRange ::  Frequency Range to filter the signal; if you do not want to filter, put []
%
% 18 juin 2015 / 13 octobre 2015




function  TriggerCSD(FileBase, fMode, State, Channels, WinSec, Name, Ref, FreqRange)

switch fMode
    
    
    case 'compOscill'
        
        eeg=TriggerCSD_load(FileBase, Channels, State, FreqRange);
        
        %%%%%%%%%%%%%%%%%
        allm = LocalMinima(eeg(:,Ref),50,0);% all local minima (allm1) avec NotCloserThan = 50 sinon ca marche pas
        figure('name','Check the Local Minima detection','NumberTitle','off')
        hold on;plot(eeg(:,Ref),'g'); plot(eeg(:,Ref));scatter(allm, ones(length(allm),1),'xm');
        Stat= input('Are Local Minima ok ? (Y=yes, N=no)  ', 's');
        if strcmp(Stat, 'N')
            disp('Change your RefCh!')
            return
        end
        clear Stat
        
        % Eliminer les outsiders et prendre les good local minima (glm1)
        figure('name','Get only the good Local Minima','NumberTitle','off')
        hist(log(abs(eeg(allm,Ref))),100)
        bound= input('Enter the Local Minima boundaries (as [x y])  ');
        glm1 = allm(log(abs(eeg(allm,Ref)))>bound(1) & log(abs(eeg(allm,Ref)))<bound(2));
        %%%%%%%%%%%%%%%%%
        
        TriggerCSD_comp(FileBase, eeg, WinSec, glm1, Name)
        
    case 'compExtTrig'
        
        eeg=TriggerCSD_load(FileBase, Channels, State, FreqRange);
        
        TriggerCSD_comp(FileBase, eeg, WinSec, Ref, Name)
        
        
        %     case 'display'        OLD WAY
        %
        %         load([FileBase '.' mfilename '.' Name '.mat'], 'Av', 'Trange');
        %
        %         figure('name', ['CSD - ' Name],'NumberTitle','off'); hold on
        %         ColorRange=PlotCSD(Av(:,:)',Trange,[],2, [], 2);
        %
        
    case 'display'
        
        load([FileBase '.' mfilename '.' Name '.mat'], 'Av', 'Trange');
        
        % from PlotCSD.m  ::
        y=Av(:,:)';
        t=Trange;
        clear Av Trange
        step=1;
        
        ch=[step+1:size(y,2)-step];
        csd = y(:,ch+step) - 2*y(:,ch) + y(:,ch-step);
        csd = -csd/(step^2);
        ColorRange = max(abs(csd(:)));
        ColorRange = [-ColorRange ColorRange];
        
        
        %  From PlotManyCh.m  ::
        x=y(:, 1+step:end-step);
        x=squeeze(x);
        nChannels = min(size(x));
        nTime = max(size(x));
        scale=2;
        if (size(x,1)~=nTime)
            x=x';
        end
        ChAmp = range(x);
        x = x*scale;
        newx = x - repmat(mean(x([1:2,end-1:end],:),1),nTime,1);
        shift = max(ChAmp);
        newx = newx - shift/2-repmat([0:nChannels-1]*shift,nTime,1);
        
        
        figure('name', [FileBase ' - CSD - ' Name],'NumberTitle','off'); hold on
        CurSrcDns(y, t, 'c', [], [], [], step, ColorRange);colorbar; xlabel('Time, ms'); ylabel('Channels')
        curAxes = get(gca,'Position');
        h = axes('Position',curAxes);
        plot(t,newx,'k', 'LineWidth', 1.5);
        set(gca, 'color', 'none');axis tight;axis off
     
 
end

function eeg=TriggerCSD_load(FileBase, Channels, State, FreqRange)

Par=LoadPar([FileBase '.xml']);
disp('Loading data...')
eeg=LoadBinary([FileBase '.eeg'], Channels, Par.nChannels);
if ~isempty(State)
    STA=load([FileBase '.sts.' State]);
    eeg=SelectPeriods(eeg(:,:), STA, 'c', 1);
end


if ~isempty(FreqRange)
    eeg = ButFilter(eeg,  2, FreqRange/(Par.lfpSampleRate/2),'bandpass');
end

function TriggerCSD_comp(FileBase, eeg, WinSec, Ref, Name)

disp('Computing CSD...')
Par=LoadPar([FileBase '.xml']);
% fenetre pour le triggered average
win=round(WinSec*Par.lfpSampleRate); % pour theta, 5000 ms

% retirer les indices de bords
a=find(Ref >win,1, 'first');
b=find(Ref +win <length(eeg), 1, 'last');
Ref=Ref(a:b);
clear a b

% trigered average
nChannels=size(eeg,2);
Av=NaN(nChannels, ((win*2)+1));

for nch=1:nChannels
    % Pour chaque canal autour des glm
    Tot=NaN(length(Ref), ((win*2)+1));
    for n=1:length(Ref)
        Tot(n,:)=eeg(Ref(n)-win : Ref(n)+win,nch) ;
    end
    Av(nch, :)=mean(Tot,1);
end

% Faire l'axe x
Srange = -win(1):win(end);
Trange = linspace(-win(1),win(end),length(Srange));
Trange=(Trange./Par.lfpSampleRate).*1000;

save([FileBase '.' mfilename '.' Name '.mat'], 'Av', 'Trange');
