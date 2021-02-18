% function  BurstIndexStates(FileBase, fMode, SubSet)
%
% Computes Burst Indices and peak burst freq values for Theta and So
%
% FileBase : Root name of the experiment in a string
% fMode: 'compute' or 'display'
% SubSet : Name it as you wish ;  'All' ; Structure Data with subfields T G THE SWS Par Map
%                !! Data.T must be in original sampling rate !!
% SpkGps : enter all the SpikeGroup numbers belonging to the Subset if chosen, ex for AnatGp2: [6:10] or [6 7 8 9 10]
% CluLocOnly : run the function only on the clusters recorded by this list of channels (from 1)

% Made according to the Map matrix you get from LoadCluRes.m
% Saves the output as FileBase.BurstIndexStates.SubSet.mat


function BurstIndexStates(FileBase, fMode, SubSet, SpkGps, CluLocOnly)



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
        
        
        
        %         BinSize = (Par.SampleRate/1000); %1 ms
        %         HalfBins = 50; % ms
        
        acgTheB=zeros(101, size(Map,1));
        acgSoB=zeros(101, size(Map,1));
        
        for n=1:size(Map,1)
            t=Tthe(Gthe==Map(n,1)); % theta
            if ~isempty(t)
                if length(t) >3
                    %             [acg, tbinB] = CCG(t, ones(size(t)), BinSize, HalfBins, Par.SampleRate,unique(ones(size(t))),'hz');
                    [acg, tbinB] = CCG(t,ones(length(t),1), 'binSize',1, 'duration', 100,'sampling', Par.SampleRate);
                    if ~isempty(acg); acgTheB(:,n)=acg; end
                end
            end
            clear t
            
            t=Tso(Gso==Map(n,1)); % so
            if ~isempty(t)
                if length(t) >3
                    %             [acg, tbinB] = CCG(t, ones(size(t)), BinSize, HalfBins, Par.SampleRate,unique(ones(size(t))),'hz');
                    [acg, tbinB] = CCG(t,ones(length(t),1), 'binSize',1, 'duration', 100,'sampling', Par.SampleRate);
                    if ~isempty(acg); acgSoB(:,n)=acg; end
                end
            end
            clear t
        end
        
        
        
        
        
        BurstsThe=zeros(4,size(Map,1)); % (1)peak, (2)baseline, (3) prep index, (4)index normalized
        BurstsSo=zeros(4,size(Map,1));
        
        for n=1:size(Map,1)
            BurstsThe(1,n)=max(acgTheB(52:61,n));
            BurstsThe(2,n)=mean(acgTheB(91:end,n));
            BurstsSo(1,n)=max(acgSoB(52:61,n));
            BurstsSo(2,n)=mean(acgSoB(91:end,n));
        end
        
        BurstsThe(3,:)=BurstsThe(1,:)-BurstsThe(2,:); %get index
        BurstsSo(3,:)=BurstsSo(1,:)-BurstsSo(2,:);
        
        x=find(BurstsSo(3,:)<0); %negative id :: id/ mean baseline
        if ~isempty(x)
            for m=1:length(x)
                BurstsSo(4,x)=BurstsSo(3,x) ./ BurstsSo(2,x);
            end
        end; clear x
        
        x=find(BurstsSo(3,:)>0); %positive id :: id/ peak
        if ~isempty(x)
            for m=1:length(x)
                BurstsSo(4,x)=BurstsSo(3,x) ./ BurstsSo(1,x);
            end
        end; clear x
        
        x=find(BurstsThe(3,:)<0); %negative id :: id/ mean baseline
        if ~isempty(x)
            for m=1:length(x)
                BurstsThe(4,x)=BurstsThe(3,x) ./ BurstsThe(2,x);
            end
        end; clear x
        
        x=find(BurstsThe(3,:)>0); %positive id :: id/ peak
        if ~isempty(x)
            for m=1:length(x)
                BurstsThe(4,x)=BurstsThe(3,x) ./ BurstsThe(1,x);
            end
        end; clear x
        
        
        if strcmp(SubSet, 'All')
            save([FileBase '.' mfilename '.' SubSet '.mat'], 'BurstsThe','BurstsSo','acgTheB','acgSoB','tbinB', 'Map');
        else
            save([FileBase '.' mfilename '.Data.mat'], 'BurstsThe','BurstsSo','acgTheB','acgSoB','tbinB', 'Map');
        end
        
        
    case 'display'
        
        
        if FileExists([FileBase '.' mfilename '.' SubSet '.mat'])
            load([FileBase '.' mfilename '.' SubSet '.mat']);
        else
            disp('Compute first !!')
        end
        
        figure('name',['Burst Index - ' SubSet],'NumberTitle','off')
        
        % Burst indices
        subplot 231; scatter(1:size(Map,1),BurstsThe(4,:),'ob');
        hold on; line([0 size(Map,1)], [0.6 0.6], 'Color', [0 1 0], 'LineStyle', '--')
        ylabel('Id Burst Hpc'); title('Theta'); ylim([-1.2 1.2]);set(gca,'XTick',[]); xlabel('neurons'); xlim([ -1 size(Map,1)+1])
        
        subplot 232; scatter(1:size(Map,1),BurstsSo(4,:),'om');
        hold on; line([0 size(Map,1)], [0.6 0.6], 'Color', [0 1 0], 'LineStyle', '--')
        title('SO'); ylim([-1.2 1.2]);set(gca,'XTick',[]);xlabel('neurons'); xlim([ -1 size(Map,1)+1])
        
        subplot 233;
        for n=1:size(Map,1)
            hold on ; scatter(1, BurstsThe(4,n),'ob');
            scatter(2, BurstsSo(4,n),'om');
            plot(1:2, [BurstsThe(4,n) BurstsSo(4,n)],'k');
        end; xlim([0.5 2.5]); ylim([-1.2 1.2]);set(gca,'XTick', [1 2]); set(gca,'XTickLabel', {'The', 'So'});
        
        % Peak Burst values
        
        xThe=find(BurstsThe(4,:) >= 0.6);
        xSo=find(BurstsThe(4,:) >= 0.6);
        
        subplot 234
        boxplot([BurstsThe(1,xThe); BurstsSo(1,xSo)]'); hold on
        
        scatter(repmat(0.7, length(xThe),1), BurstsThe(1,xThe),'MarkerEdgeColor',[0 0 1]); % 'MarkerEdgeColor',[0.8,0.2,0.1]
        scatter(repmat(1.7, length(xSo),1), BurstsSo(1,xSo),'MarkerEdgeColor',[1 0 1]); % 'MarkerEdgeColor',[0.1,0.8,0.2]
        ylabel('Burst Peak Frequency for Bursting Neurons Only, Hz')
        set(gca,'XTick',1:2);
        set(gca,'XTickLabel', {'THE', 'SO'});
        %
        %         subplot 235
        %         for n=1:size(BurstsThe,2)
        %             hold on
        %             scatter(1, BurstsThe(1,n),'MarkerEdgeColor',[0.8,0.2,0.1],'MarkerFaceColor',[0.8,0.2,0.1]);
        %             scatter(2, BurstsSo(1,n),'MarkerEdgeColor',[0.1,0.8,0.2],'MarkerFaceColor',[0.1,0.8,0.2]);
        %             plot(1:2, [BurstsThe(1,n); BurstsSo(1,n)],'Color',[0.2,0.2,0.2])
        %
        %             xlim([0.5 2.5]);
        %             ylabel('Burst Peak Frequency, Hz')
        %             set(gca,'XTick',1:2);
        %             set(gca,'XTickLabel', {'THE', 'SO'});
        %         end
        %
        
end


