% function Comodulogram(FileBase, fMode, State, Channels, Sets)
% Sets = example 'mPfcHpc';
%


function Comodulogram(FileBase, fMode, State, Channels, Sets)


switch fMode
    case 'compute'
        
        Par=LoadPar([FileBase '.xml']);
        STA = load([FileBase '.sts.' State]);
        eeg = LoadBinary([FileBase '.eeg'], Channels, Par.nChannels);
        eeg = SelectPeriods(eeg(:,:),STA,'c',1);
        
        
        WinLengthSec=4;
        WinLengthSample = 2^round(log2(WinLengthSec*Par.lfpSampleRate));
        nFFT = 2*WinLengthSample;
        nChannels=length(Channels);
        
        if nChannels ==2
            
            [co,CoM.f, Pval, CoRlo, CoRup] = Comodugram(eeg(:,:), nFFT, Par.lfpSampleRate, [1 150], WinLengthSample, 3,'linear', 1);
            CoM.co(:,:)=co(:,:,1,2); clear co
            CoM.Pval(:,:)=Pval(:,:,1,2); clear co
            CoM.CoRlo(:,:)=CoRlo(:,:,1,2); clear co
            CoM.CoRup(:,:)=CoRup(:,:,1,2); clear co
            
        else
            for n=1:nChannels-1
                for nn=n+1:nChannels
                    [co,CoM.f, Pval, CoRlo, CoRup] = Comodugram(eeg(:,[n nn]), nFFT, Par.lfpSampleRate, [1 150], WinLengthSample, 3,'linear', 1);
                    CoM.co(:,:,n,nn)=co(:,:,n,nn); clear co
                    CoM.Pval(:,:,n,nn)=Pval(:,:,n,nn); clear co
                    CoM.CoRlo(:,:,n,nn)=CoRlo(:,:,n,nn); clear co
                    CoM.CoRup(:,:,n,nn)=CoRup(:,:,n,nn); clear co
                end
            end
        end
        
        
        save([FileBase '.' mfilename '.' Sets '.' State '.mat'], 'CoM');
        
    case 'display'
        
        load([FileBase '.' mfilename '.' Sets '.' State '.mat']);
        
        figure('name', [FileBase 'Commodugram ' State],'NumberTitle','off')
        imagesc(CoM.f, CoM.f, CoM.co); axis xy;
        colormap('Jet'); colorbar
        figure('name', [FileBase 'Commodugram ' State ' - pvalues'],'NumberTitle','off')
        imagesc(f,f,Pval(:,:,1,2));set(gca,'ydir','norm'); coloration('pvalue');
        
        
    case 'displaytout'
        
        load([FileBase '.' mfilename '.' Sets '.' State '.mat']);
        
        figure('name', [FileBase 'Commodugram ' State],'NumberTitle','off')
        for n=1:length(Channels)
            subplotfit(n, length(Channels))
            imagesc(CoM.f, CoM.f, CoM.co(:,:,n)); axis xy;
            %         xlabel(['Channel ' num2str(Channels(n))]); ylabel(['Channel ' num2str(Channels(3))]) ;
            colormap('Jet'); colorbar
        end
        
    case 'displaytoutfig'
        
        load([FileBase '.' mfilename '.' Sets '.' State '.mat']);
        
        figure('name', [FileBase 'Commodugram ' State],'NumberTitle','off')
        imagesc(CoM.f, CoM.f, CoM.co(:,:,1)); axis xy;
        xlabel(['mPFC']); ylabel(['Hpc']) ; colormap('Jet'); colorbar
        clim([-0.3 0.3]);
        
end


