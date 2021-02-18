% function OutArgs = Gamma2Units(FileBase, fMode, WhereGam, SubSet, SpkGps, CluLocOnly, State, Normalization)
%
% FileBase : Root name of the experiment in a string
% fMode: 'compute' or 'display'
% WhereGam :: where gamma has been detected (as in DetectGammaBursts.m)
% SubSet : Name it as you wish or 'All' ;
% SpkGps : enter all the SpikeGroup numbers belonging to the Subset if chosen, ex for AnatGp2: [6:10] or [6 7 8 9 10]
% CluLocOnly : run the function only on the clusters recorded by this list of channels (from 1)
% State : 'THE, 'SWS'....
%
% updated July 5th 2016 with no normalization (because GBtime_new is
% beginning of the burst) and GBtime_new (all good gamma events) 
% AND jitter corrected!

function OutArgs = Gamma2Units(FileBase, fMode, WhereGam, SubSet, SpkGps, CluLocOnly, State)

switch fMode
    case 'compute'
        Par = LoadPar([FileBase '.xml']);

        load([FileBase '.DetectGammaBursts.' WhereGam '.' State '.mat'], 'GBtime_new'); % for a given state
        
        %~% Loading units %~%
        
        if strcmp(SubSet, 'All');
            if ~FileExists([FileBase '.CluRes.mat'])
                [T,G,Map]=LoadCluRes(FileBase);
            else
                load([FileBase '.CluRes.mat']);
            end
        else
            [T,G,Map]=LoadCluRes(FileBase, SpkGps,0, CluLocOnly);
        end
        
        T =  max(1,round(T*Par.lfpSampleRate/Par.SampleRate)); % get spike times to the sample rate of phase.avoid zero
        
        
        %~% Get good spikes in the given state %~%
 
        if ~exist([FileBase '.sts.' State],'file')
            disp('No sts file for this state');return
        end
        STA = load([FileBase '.sts.' State]);
        
        % [T, Ind] = SelPerDiscr(T, STA, 1, 1); Cette fonction: c'est de la merde!! pas de respect des timestamps originaux
        [T, Ind] = SelPerDiscr2(T, STA);
        G = G(Ind);
        OutArgs.map=Map;
        clear Ind
    
        %~% Prep for Firing proba
        
        if FileExists([FileBase '.FiringRateStates.' SubSet '.mat'])
            load([FileBase '.FiringRateStates.' SubSet '.mat'], 'FRate');
        else
            warning(['You must run FiringRateStates.m on ' SubSet ' first!']);
            return
        end
        OutArgs.FProb=[];

        
        %~% CCG and significance
        
        binsize=10; halbins=50;
        [ccg, OutArgs.tbin] = Trains2CCG({GBtime_new,T},{1,G},binsize,halbins,Par.lfpSampleRate,'count');
        OutArgs.Ccg = squeeze(ccg(:,1,2:end));
        
        OutArgs.FProb=NaN(size(OutArgs.Ccg));
        for n=1:size(OutArgs.FProb,2)
            OutArgs.FProb(:,n) =  OutArgs.Ccg(:,n)./FRate(2,n);
        end
        
        
        
                %  CCG for jittering data
        njitter=1000;jscale=100; alpha=0.05;
        disp('... computing significance, may take time ...'); tic
         
        for n=1:size(OutArgs.map,1)
            
%             [ccg, tbin] = CCG([GBtime_new;T(G==n)],[ones(size(GBtime_new));2*ones(size(T(G==n)))], binsize, ...
%                 halbins, Par.lfpSampleRate,[1,2],'count');
%             OutArgs.Ccg(:,n)=squeeze(ccg(:,1,2));
%             OutArgs.tbin=tbin;
            disp(['-> computing neuron ' num2str(n) '/' num2str(size(OutArgs.map,1))])
%             
%                 % Firing Proba
%             OutArgs.FProb= [OutArgs.FProb OutArgs.Ccg(:,n)./FRate(2,n)];
            
                %  CCG for jittering data
            for x=1:njitter
                T_jitter = T(G==n) + round(2*round(Par.lfpSampleRate/1000*jscale)*rand(size(T(G==n)))-1 ...
                    *round(Par.lfpSampleRate/1000)*jscale);
                [jccg, tJ] = CCG([GBtime_new;T_jitter],[ones(size(GBtime_new));2*ones(size(T_jitter))], binsize, ...
                    halbins, Par.lfpSampleRate,[1,2],'count');
                ccgj(:,x)=jccg(:,1,2);
            end
            
                %  Computes the pointwise line
            ccgjm(:,n)  = mean(ccgj,2);
            ccgjmax(:,n)=ccgjm(n) + 2*std(ccgj,[],2);
            ccgjmin(:,n)=ccgjm(n) - 2*std(ccgj,[],2);
 
%                 signifpoint = njitter*alpha;
%                 for n=1:length(tJ)
%                   sortjitterDescend  = sort(jccg(n,:),'descend');
%                   sortjitterAscend   = sort(jccg(n,:),'ascend');
%                   ccgjmax(n) = sortjitterDescend(signifpoint);
%                   ccgjmin(n) = sortjitterAscend(signifpoint);
%                 end

    
    
            
            OutArgs.ccgjm=ccgjm;
            OutArgs.ccgjptMax=ccgjmax;
            OutArgs.ccgjptMin=ccgjmin;
            OutArgs.njitter=njitter;
            OutArgs.jscale=jscale;
            
        end
        

        save([FileBase '.' mfilename '.' WhereGam '.' SubSet  '.' State '.mat'], 'OutArgs');
        
    case 'display'
        load([FileBase '.' mfilename '.' WhereGam '.' SubSet  '.' State '.mat'], 'OutArgs');
        
        if ~isempty(OutArgs) %#ok<*NODEF>
            figure('name',['Gamma2Units - ' WhereGam ' - ' SubSet ' - ' State],'NumberTitle','off')
            set(0, 'DefaultFigurePosition', [-1118         483        1022         776])
            ncells =size(OutArgs.Ccg,2);
            for n=1:ncells
                subplotfit(n,ncells);
                
                bar(OutArgs.tbin, squeeze(OutArgs.Ccg(:,n))'); hold on; axis tight
                line(OutArgs.tbin,OutArgs.ccgjm(:,n),'linestyle','--','color','b','LineWidth',1.5)
                line(OutArgs.tbin,OutArgs.ccgjptMax(:,n),'linestyle','--','color','r','LineWidth',1.5)
                line(OutArgs.tbin,OutArgs.ccgjptMin(:,n),'linestyle','--','color','r','LineWidth',1.5)
                title(num2str(OutArgs.map(n,2:3)));
            end
        else
            return
        end
        
    case 'dispFProb'
        
        load([FileBase '.' mfilename '.' WhereGam '.' SubSet  '.' State '.mat']);
        if ~isempty(OutArgs)
            figure('name',['Gamma2Units - ' WhereGam ' - ' SubSet ' - ' State],'NumberTitle','off')
            ncells =size(OutArgs.ccg,2);
            for n=1:ncells
                subplotfit(n,ncells);
                bar(OutArgs.tbin, squeeze(OutArgs.FProb(:,n))');
                line([0 0], [0, max(OutArgs.FProb(:,n))], 'Color', [1 0 0], 'LineWidth', 2)
                axis tight
                title(num2str(OutArgs.map(n,2:3)));
            end
        else
            return
        end
        
end

%% OLD VERSION

% function OutArgs = Gamma2Units(FileBase, fMode, WhereGam, SubSet, SpkGps, CluLocOnly, State)
% 
% switch fMode
%     case 'compute'
%         Par = LoadPar([FileBase '.xml']);
%         %         if nargin < 8
%         %             Normalization = 0;
%         %         end
%         
%         %         load([FileBase '.DetectGammaBursts.' WhereGam '.' State '.mat'], 'GBtime_new', 'GBduration'); % for a given state
%         load([FileBase '.DetectGammaBursts.' WhereGam '.' State '.mat'], 'GBtime_new'); % for a given state
%         
%         % % Loading units
%         
%         if strcmp(SubSet, 'All');
%             if ~FileExists([FileBase '.CluRes.mat'])
%                 [T,G,Map]=LoadCluRes(FileBase);
%             else
%                 load([FileBase '.CluRes.mat']);
%             end
%         else
%             [T,G,Map]=LoadCluRes(FileBase, SpkGps,0, CluLocOnly);
%         end
%         
%         T =  max(1,round(T*Par.lfpSampleRate/Par.SampleRate)); % get spike times to the sample rate of phase.avoid zero
%         
%         
%         if ~exist([FileBase '.sts.' State],'file')
%             disp('No sts file for this state');return
%         end
%         STA = load([FileBase '.sts.' State]);
%         
%         %         [Ts, Ind] = SelPerDiscr(T, STA, 1, 1); Cette fonction: c'est de la merde!! pas de respect des timestamps originaux
%         [T, Ind] = SelPerDiscr2(T, STA);
%         G = G(Ind);
%         clear Ind
%         
%         
%         %
%         %         if Normalization == 0
%         
%         [ccg, OutArgs.tbin] = Trains2CCG({GBtime_new,T},{1,G},10,50,Par.lfpSampleRate,'scale');
%         OutArgs.ccg = squeeze(ccg(:,1,2:end));
%         OutArgs.map = [Map NaN(size(Map,1),2)];
%         OutArgs.norm=0;
%         
% 
% 
%         if FileExists([FileBase '.FiringRateStates.' SubSet '.mat'])
%             load([FileBase '.FiringRateStates.' SubSet '.mat'], 'FRate');
%         else
%             warning(['You must run FiringRateStates.m on ' SubSet ' first!']);
%             return
%         end
%         
%         OutArgs.FProb=NaN(size(OutArgs.ccg));
%         for n=1:size(OutArgs.FProb,2)
%             OutArgs.FProb(:,n) =  OutArgs.ccg(:,n)./FRate(2,n);
%         end
%         
%         %%% Compute significance
%         
%         %  CCG for jittering data
%         njitter=500;jscale=100;
%         disp('... computing significance, may take time ...'); tic
%         
%         
%         for j=1:njitter
%             T_jitter = T + 2*(20*jscale)*rand(size(T))-1*20*jscale;
%             ccg = CCG([GBtime_new;T_jitter],[ones(size(GBtime_new));2*ones(size(T_jitter))], 10, 50, Par.lfpSampleRate,[1,2],'scale');
%             ccgj(:,j)=ccg(:,1,2); %#ok<*AGROW>
%             ccgjmax(j)=max(ccgj(:,j));
%             ccgjmin(j)=min(ccgj(:,j));
%         end
%         clear ccg
%         OutArgs.ccgjm  = mean(ccgj,2); toc
%         
%         %  Computes the pointwise line
%         alph=0.05;
%         signifpoint = njitter*alph;
%         for j=1:length(OutArgs.tbin)
%             sortjitterDescend  = sort(ccgj(j,:),'descend');
%             sortjitterAscend   = sort(ccgj(j,:),'ascend');
%             ccgjptMax(j) = sortjitterDescend(signifpoint);
%             ccgjptMin(j) = sortjitterAscend(signifpoint);
%         end
%         
%         OutArgs.ccgjptMin=ccgjptMin;
%         OutArgs.ccgjptMax=ccgjptMax;
%         clear njitter jscale alph sortjitterDescend sortjitterAscend ccgjptMax ccgjptMin ccgjm ccg ccgjmax ccgjmin ccgj  t
%         
%         %
%         %         else
%         %
%         %             uclu=unique(G);
%         %             OutArgs.ccg=[];
%         %             OutArgs.tbin=[];
%         %
%         %             for clu=1:length(uclu) % for each neuron
%         %                 t=T(G==uclu(clu));
%         %                 [tsOffsets, ts1idx, ts2idx] = crosscorrelogram(GBtime_new/Par.lfpSampleRate, t/Par.lfpSampleRate, [-0.5 0.5]);
%         %                  tsOffsetsN=NaN(length(tsOffsets),1);
%         %                 for n=1:length(tsOffsets)
%         %                     tsOffsetsN(n)= tsOffsets(n)/GBduration(ts1idx(n));
%         %                 end
%         %                 [ccg, tbin]=hist(tsOffsetsN, 200);
%         %                 OutArgs.ccg=[OutArgs.ccg ; ccg];
%         %                 OutArgs.Offsets{clu}= tsOffsetsN;
%         %                 OutArgs.tbin=[OutArgs.tbin ; tbin];
%         %             end
%         %
%         %             OutArgs.ccg=OutArgs.ccg';
%         %             OutArgs.norm=1;
%         %             OutArgs.map = Map;
%         %
%         
%         %         end
%         
%         save([FileBase '.' mfilename '.' WhereGam '.' SubSet  '.' State '.mat'], 'OutArgs');
%         
%     case 'display'
%         load([FileBase '.' mfilename '.' WhereGam '.' SubSet  '.' State '.mat'], 'OutArgs');
%         
%         if ~isempty(OutArgs) %#ok<*NODEF>
%             figure('name',['Gamma2Units - ' WhereGam ' - ' SubSet ' - ' State],'NumberTitle','off')
%             set(0, 'DefaultFigurePosition', [-1118         483        1022         776])
%             ncells =size(OutArgs.ccg,2);
%             for ii=1:ncells
%                 subplotfit(ii,ncells);
%                 %                 bar(OutArgs.tbin(ii,:), squeeze(OutArgs.ccg(:,ii))');
%                 bar(OutArgs.tbin, squeeze(OutArgs.ccg(:,ii))'); hold on; axis tight
%                 line([0 0], [0, max(OutArgs.ccg(:,ii))], 'Color', [1 0 0], 'LineWidth', 2)
%                 line([-1 -1], [0 max(OutArgs.ccg(:,ii))],'linestyle','--','Color',[1 0.5 0],'LineWidth',0.5)
%                 line([1 1], [0 max(OutArgs.ccg(:,ii))],'linestyle','--','Color',[1 0.5 0],'LineWidth',0.5)
%                 %                  xlim([-5 5])
%                 title(num2str(OutArgs.map(ii,2:3)));
%             end
%         else
%             return
%         end
%         
%     case 'dispFProb'
%         
%         load([FileBase '.' mfilename '.' WhereGam '.' SubSet  '.' State '.mat']);
%         if ~isempty(OutArgs)
%             figure('name',['Gamma2Units - ' WhereGam ' - ' SubSet ' - ' State],'NumberTitle','off')
%             ncells =size(OutArgs.ccg,2);
%             for ii=1:ncells
%                 subplotfit(ii,ncells);
%                 bar(OutArgs.tbin, squeeze(OutArgs.FProb(:,ii))');
%                 line([0 0], [0, max(OutArgs.FProb(:,ii))], 'Color', [1 0 0], 'LineWidth', 2)
%                 axis tight
%                 title(num2str(OutArgs.map(ii,2:3)));
%             end
%         else
%             return
%         end
%         
% end

