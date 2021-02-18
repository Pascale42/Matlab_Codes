% myCSD(FileBase, eeg, fMode, What, refch, FreqRange, SR);
%
% eeg :: can be the matrix to compute CSD  ! OR ! can be the extension of a
%           file to load; 
%           ex: FileBase.lfp.THE-HPC.mat, here eeg= 'lfp.THE-HPC'
% fMode :: 'compute' or 'display'
% What :: extension name for saving the output
%           ex: What='THE-Hpc' will give FileBase.CSD.THE-Hpc.mat
% refch :: reference channel for the trigger average (largest power of your signal of interest)
% FreqRange :: [x xx];
% RS :: Sampling Rate
%

function myCSD(FileBase, eeg, fMode, What, refch, FreqRange, SR)

switch fMode
    
    case 'compute'
        
        if ischar(eeg)
            load([FileBase '.' eeg '.mat'], 'eeg');
        end
        
        %Filtrer
        disp('... filtering...')
        feeg = ButFilter(eeg, 2, FreqRange/(SR/2),'bandpass');
        
        
        % all local minima (allm1) avec NotCloserThan = 50 sinon ca marche pas
        allm = LocalMinima(feeg(:,refch),50,0);
        
        %figure; hold on;plot(eeg(:,refch),'g'); plot(fvHpcT(:,refch));scatter(allm, ones(length(allm),1),'xm');
        clear eeg
        
        %Eliminer les outsiders et prendre les good local minima (glm1)
        figure('name', 'Get the boundaries of the local minima distribution','NumberTitle','off');
        hist(log(abs(feeg(allm,refch))),100)
        [x,y]=ginput(2); 
        close('name', 'Get the boundaries of the local minima distribution','NumberTitle','off');
        glm1 = allm(log(abs(feeg(allm,refch)))>x(1) & log(abs(feeg(allm,refch)))<x(2));
        clear xy
        
        %fenetre pour le triggered average
        f=(FreqRange(1)+FreqRange(2))/2;
        f=round(1000/f);
        win=round(f/1000*SR); % pour theta, 250 ms
        clear f
        
        % retirer les indices de bords
        a=find(glm1 >win,1, 'first');
        b=find(glm1 +win <length(feeg), 1, 'last');
        glm=glm1(a:b);
        clear a b glm1
        
        
        % trigered average
        nChannels=size(feeg,2);
        Av=NaN(nChannels, ((win*2)+1));
        
        %csd = CSD(eegrHT);
         disp('... computes trigered average...')
        for ch=1:nChannels
            % Pour chaque canal, matrice avec lfp cut autour des glm
            Tot=NaN(length(glm), ((win*2)+1));
            for n=1:length(glm)
                Tot(n,:)=feeg([glm(n)-win : glm(n)+win],ch) ;
            end
            Av(ch, :)=mean(Tot,1);
        end
        clear Tot glm
        
        % Faire l'axe x
        Srange = [-win(1):win(end)];
        Trange = linspace(-win(1),win(end),length(Srange));
        Trange=Trange*2;
        clear Srange
        
        save([FileBase '.CSD.' What '.mat'], 'Av', 'Trange', 'win',  'feeg');
        
        
    case 'display'
        
        load([FileBase '.CSD.' What '.mat'], 'Av', 'Trange');
        
        figure('name', ['CSD ' What],'NumberTitle','off');
        PlotCSD(Av(:,:)',Trange,[],2, [], 2); colormap('jet')
        
end
