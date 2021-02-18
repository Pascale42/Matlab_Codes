% function FRate = FiringRateStates(FileBase, fMode, SubSet, SpkGps, CluLocOnly)
%
% Computes Mean Firing Rate for Theta and So, and overall
%
% FileBase : Root name of the experiment in a string
% fMode: 'compute' or  'display'
% SubSet : Name it as you wish ;  'All' ; Structure Data with subfields T G THE SWS Par Map
%                !! Data.T must be in original sampling rate !!
% SpkGps : enter all the SpikeGroup numbers belonging to the Subset if chosen, ex for AnatGp2: [6:10] or [6 7 8 9 10]
% CluLocOnly : run the function only on the clusters recorded by this list of channels (from 1)
%
% The output is FRate:
% It gives for each column (i.e. neuron) the mean firing rate
% for THE (first line) and for SWS (second line)
% Last line is the global (no brainstate) mean firing rate
% Made according to the Map matrix you get from LoadCluRes.m
% Saves the output as FileBase.FiringRateStates.SubSet.mat


function FRate = FiringRateStates(FileBase, fMode, SubSet, SpkGps, CluLocOnly)



switch fMode
    
    case 'compute'
        
        % Loading units
       
        if ischar(SubSet)  %%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if strcmp(SubSet, 'All')
                if ~FileExists([FileBase '.CluRes.mat'])
                    [T,G,Map, Par]=LoadCluRes(FileBase);
                else
                    load([FileBase '.CluRes.mat']);
                end
            else
                [T,G,Map,Par]=LoadCluRes(FileBase, SpkGps,0, CluLocOnly);
            end
                        
            % Theta
            if FileExists([FileBase '.sts.THE'])
                THE=load([FileBase '.sts.THE']);
                THE=round(THE*(Par.SampleRate/Par.lfpSampleRate)); % put the same sampling rate as spikes
                [Tthe, Indthe] = SelPerDiscr(T, THE, 1, 1);
                Gthe = G(Indthe);
            else
                Gthe=[]; Tthe=[];
            end
            
            % SO
            if FileExists([FileBase '.sts.SWS'])
                SO=load([FileBase '.sts.SWS']);
                SO=round(SO*(Par.SampleRate/Par.lfpSampleRate));
                [Tso, Indso] = SelPerDiscr(T, SO, 1, 1);
                Gso = G(Indso);
            else
                Gso=[]; Tso=[];
            end
            
        else %%%%%%%%%%%%%%%%%%%%%%%%%%
            
            T=SubSet.T;
            G=SubSet.G;
            Map=SubSet.Map;
            Par=SubSet.Par;
            
            if  isfield(SubSet, 'SWS')
                SO=SubSet.SWS;
                SO=round(SO*(Par.SampleRate/Par.lfpSampleRate));
                [Tso, Indso] = SelPerDiscr(T, SO, 1, 1);
                Gso = G(Indso);
            else
                Gso=[]; Tso=[];
            end
            
            if  isfield(SubSet, 'THE')
                THE=SubSet.THE;
                THE=round(THE*(Par.SampleRate/Par.lfpSampleRate)); 
                [Tthe, Indthe] = SelPerDiscr(T, THE, 1, 1);
                Gthe = G(Indthe);
            else
                Gthe=[]; Tthe=[];
            end
            
            
        end %%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        
        
        MinISI = 5; % in samples; avoids double detected spikes if any
        FRate=zeros(3,size(Map,1)); % for all neurons (rows), in The (line1), So (line2) and All (line3)
        
        for n=1:size(Map,1)
            % THE
            myClu = Gthe==Map(n,1);
            myRes = Tthe(myClu);
            myISI = diff(myRes);
            goodISI = myISI>MinISI;
            FRate(1,n) = Par.SampleRate / mean(myISI(goodISI));
            % SWS
            myClu = Gso==Map(n,1);
            myRes = Tso(myClu);
            myISI = diff(myRes);
            goodISI = myISI>MinISI;
            FRate(2,n) = Par.SampleRate / mean(myISI(goodISI));
            % ALL
            myClu = G==Map(n,1);
            myRes = T(myClu);
            myISI = diff(myRes);
            goodISI = myISI>MinISI;
            FRate(3,n) = Par.SampleRate / mean(myISI(goodISI));
        end
        
        
        if strcmp(SubSet, 'All')
            save([FileBase '.' mfilename '.' SubSet '.mat'], 'FRate', 'Map');
        else
            save([FileBase '.' mfilename '.Data.mat'], 'FRate', 'Map');
        end
        
      
        
        
    case 'display'
        
        
        if FileExists([FileBase '.' mfilename '.' SubSet '.mat'])
            load([FileBase '.' mfilename '.' SubSet '.mat']);
        else
            disp('Compute first !!')
        end
        
        figure('name',['Mean Firing Rates - ' SubSet],'NumberTitle','off')
        
        subplot(1,3,1:2)
        % plots statistic boxes
        boxplot(FRate(:,:)','colors',[0.8,0.2,0.1;0.1,0.8,0.2;0.2,0.2,0.2],'notch','on'); hold on
        % plot individual values
        scatter(repmat(0.7, size(FRate,2),1), FRate(1,:),'MarkerEdgeColor',[0.8,0.2,0.1]);
        scatter(repmat(1.7, size(FRate,2),1), FRate(2,:),'MarkerEdgeColor',[0.1,0.8,0.2]);
        scatter(repmat(2.7, size(FRate,2),1), FRate(3,:),'MarkerEdgeColor',[0.2,0.2,0.2]);
        ylabel('Firing Rate, Hz')
        set(gca,'XTick',1:3);
        set(gca,'XTickLabel', {'THE', 'SO','ALL'});
        
        subplot(1,3,3)
        for n=1:size(FRate,2)
            hold on
            scatter(1, FRate(1,n),'MarkerEdgeColor',[0.8,0.2,0.1],'MarkerFaceColor',[0.8,0.2,0.1]);
            scatter(2, FRate(2,n),'MarkerEdgeColor',[0.1,0.8,0.2],'MarkerFaceColor',[0.1,0.8,0.2]);
            plot(1:2, FRate(1:2,n),'Color',[0.2,0.2,0.2])
            xlim([0.5 2.5]);
            ylabel('Firing Rate, Hz')
            set(gca,'XTick',1:2);
            set(gca,'XTickLabel', {'THE', 'SO'});
        end
        
end

