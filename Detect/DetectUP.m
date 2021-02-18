% function DetectUP(FileBase, Mode, Threshold)
%
% :::::::::: INPUTS ::::::::::
% - FileBase :: Root name of the file (like the eeg file)
% - Mode :: 3 possible modes:
%                - 'compute' 
%                - 'showdetection'  to display the data with the detected UPstates (green cross=beginning; pink cross=end)
%                - 'display'  to display and get the distribution of UPstates duration
% - Threshold :: optional input, fixed at 0.15 if not given
%
% :::::::::: OUTPUTS ::::::::::
%

function DetectUP(FileBase, Mode, Threshold) %#ok<INUSD>

if nargin<3
    Threshold=0.15;
end


switch Mode
    case 'compute'
        
        
        % % Loading data & Checking that all files are ready
         
        % Checks if the eeg file exists
        if FileExists([FileBase '.eeg']);
            data=LoadBinary([FileBase '.eeg'], 1, 1);
        else
            disp('Please use MakeEeg.m First!'); return
        end      
        
        % Checks if the states files exist
        States={'CTR'; 'ODN'};
        if ~FileExists([FileBase '.sts.' States{1}]) || ~FileExists([FileBase '.sts.' States{2}])
            disp('Please use MakeStates.m First!'); return
        end
        
        % Checks if you want to overwrite if UPstates already detected
        if FileExists([FileBase '.' mfilename '.mat']);
            Stat = input('Do you want to recompute this? (Y=yes, N=no)  ', 's');
            if strcmp(Stat, 'N')
                disp('UPstates are NOT redetected!'); return
            elseif strcmp(Stat, 'Y')
                disp('UPstates being redetected...');
            else
                return
            end
        end
        
        Par = LoadPar([FileBase '.xml']);
        
        
        % % Detect UP states for CTR and ODN states
        
        for s=1:length(States)
            
            % Assign data
            STA=load([FileBase '.sts.' States{s}]);
            eeg=SelectPeriods(data(:,:), STA, 'c', 1);
            
            % Get Standard deviation and offset data
            Std=std(eeg);
            eegoff=eeg-Std;
            
            % Filter eeg
            feegoff= ButFilter(eegoff, 2, [0.5 2]/(Par.lfpSampleRate/2),'bandpass');
            Stdoff=std(feegoff);
            Threshold=0.25;
            
            % SchmittTrigger
            STthresh=-Stdoff*Threshold;
            n = length(feegoff);
            
            PrevVal = [(STthresh+STthresh)/2; feegoff(1:n-1)];
            
            % Beginning/End of UPstates
            UpCrossings = find(PrevVal<STthresh & feegoff>=STthresh);
            DownCrossings = find(PrevVal>STthresh & feegoff<=STthresh);
            
            % OUTPUTS
            UPstates.int.(States{s}) = [DownCrossings UpCrossings]; % beginning/end of eachUPstates
            UPstates.dur.(States{s}) = (UpCrossings - DownCrossings)/Par.lfpSampleRate *1000; % Duration in msec
        end
        
        save([FileBase '.' mfilename '.mat'], 'UPstates');
        
    case 'showdetection'
        
        if FileExists([FileBase '.' mfilename '.mat']);
            load([FileBase '.' mfilename '.mat']);
        else
            disp('Please Compute First!'); return
        end
        
        data=LoadBinary([FileBase '.eeg'], 1, 1);
        CTR=load([FileBase '.sts.CTR']);
        eegCTR=SelectPeriods(data(:,:), CTR, 'c', 1);
        ODN=load([FileBase '.sts.ODN']);
        eegODN=SelectPeriods(data(:,:), ODN, 'c', 1);
        
        figure('name','UP States detection','NumberTitle','off')
        subplot 211; plot(eegCTR, 'k'); hold on; axis tight; title('eeg CTR'); xlabel('Time')
        scatter(UPstates.int.CTR(:,1), repmat(3e03,  length(UPstates.int.CTR(:,1)), 1), 'xg', 'LineWidth', 1.2)
        scatter(UPstates.int.CTR(:,2), repmat(3e03, length(UPstates.int.CTR(:,2)), 1), 'xm', 'LineWidth', 1.2)
        subplot 212; plot(eegODN, 'k'); hold on; axis tight; title('eeg ODN'); xlabel('Time')
        scatter(UPstates.int.ODN(:,1), repmat(3e03,  length(UPstates.int.ODN(:,1)), 1), 'xg', 'LineWidth', 1.2)
        scatter(UPstates.int.ODN(:,2), repmat(3e03, length(UPstates.int.ODN(:,2)), 1), 'xm', 'LineWidth', 1.2)
        
    case 'display'
        
        if FileExists([FileBase '.' mfilename '.mat']);
            load([FileBase '.' mfilename '.mat']);
        else
            disp('Please Compute First!'); return
        end
        
        
        figure('name','UP States duration distributions','NumberTitle','off'); % MAKE SIMILAR BINS FOR HISTOGRAMS !!!!!!
        subplot 131; Hctr=histfit(UPstates.dur.CTR,100,'kernel'); title('CTR')
        subplot 132; Hodn=histfit(UPstates.dur.ODN,100,'kernel'); title('ODN')
        
        A=who('-file',[FileBase '.' mfilename '.mat']);
        if ~sum(strcmp(A, 'FitCtr'));  
        % Get the histogram datapoints
        FitCtr(:,1)=get(Hctr(2),'XData');
        FitCtr(:,2)=get(Hctr(2),'YData');
        FitOdn(:,1)=get(Hodn(2),'XData');
        FitOdn(:,2)=get(Hodn(2),'YData');

        x=get(Hctr(1),'XData'); HistCtr(:,1)=x(1,:);
        y=get(Hctr(1),'YData'); HistCtr(:,2)=y(2,:);
        x=get(Hodn(1),'XData'); HistOdn(:,1)=x(1,:);
        y=get(Hodn(1),'YData'); HistOdn(:,2)=y(2,:);
        
        
       [MaxCtr,id]=max(FitCtr(:,2)/length(UPstates.dur.CTR));
       MaxCtr(2)=FitCtr(id,1);
       [MaxOdn,id]=max(FitOdn(:,2)/length(UPstates.dur.ODN));
       MaxOdn(2)=FitOdn(id,1); 
       
       % Kolmogorov-Smirnov test;
       % h=1; les 2 distrib sont differentes
       [h, p] = kstest2(UPstates.dur.CTR, UPstates.dur.ODN);
       
        save([FileBase '.' mfilename '.mat'], 'FitCtr', 'FitOdn', 'HistCtr', 'HistOdn', 'h', 'p', '-append');  
        end
        
        subplot 133; plot(FitCtr(:,1), FitCtr(:,2)/length(UPstates.dur.CTR),'k', 'LineWidth', 1.5); hold on; 
        plot(FitOdn(:,1), FitOdn(:,2)/length(UPstates.dur.ODN),'r', 'LineWidth', 1.5); 
        legend('CTR', 'ODN'); xlabel('Time, ms'); ylabel('Normalized count');
        if h
            title(['Distributions Normalized :: KStest = different distrib; p=' num2str(p)]);
        else
           title(['Distributions Normalized :: KStest = no difference; p=' num2str(p)]);
        end
            
        disp(['Median UPstates CTR duration = ' num2str(median(UPstates.dur.CTR)) ' ms | min ' num2str(min(UPstates.dur.CTR)) ', max ' num2str(max(UPstates.dur.CTR))])
        disp(['Median UPstates ODN duration = ' num2str(median(UPstates.dur.ODN)) ' ms | min ' num2str(min(UPstates.dur.ODN)) ', max ' num2str(max(UPstates.dur.ODN))])
end
