% function ACGsStates(FileBase, fMode, HalfBins, SubSet)
%
% Computes and Display ACGs for THE vs SWS
%
% FileBase : Root name of the experiment in a string
% fMode: 'compute' or 'display' (click on the mouse or keyboard to display the next neuron)
% HalfBins : in msec (like in Klusters)
% SubSet : 'AnatGp1', 'AnatGp2', 'All'  corresponds to the spikes of the
% first AnatGroup (Hippocampus), or of the second, or all of them
%
% Made according to the Map matrix you get from LoadCluRes.m
% Saves the output as FileBase.ACGsStates.SubSet.HalfBins.mat



function ACGsStates(FileBase, fMode, HalfBins, SubSet)



switch fMode
    case 'compute'
        
        
        if FileExists([FileBase '.CluRes.mat']);
            load([FileBase '.CluRes.mat']);
        else
            LoadCluRes(FileBase);
            load([FileBase '.CluRes.mat']);
        end
        
        disp(' ... Computing...')
        
        % Theta
        THE=load([FileBase '.sts.THE']);
        THE=round(THE*(Par.SampleRate/Par.lfpSampleRate)); % put the same sampling rate like spikes
        [Tthe Indthe] = SelPerDiscr(T, THE, 1, 1);
        Gthe = G(Indthe);
        
        % SO
        SO=load([FileBase '.sts.SWS']);
        SO=round(SO*(Par.SampleRate/Par.lfpSampleRate));
        [Tso Indso] = SelPerDiscr(T, SO, 1, 1);
        Gso = G(Indso);
        
        switch SubSet
            
            case 'AnatGp1'
                
                clu=find(Map(:,2) <= 5, 1,'last');
                Map=Map(1:clu, :);
                
                % Theta
                indx=find(Gthe<=clu);
                Gthe=Gthe(indx);
                Tthe=Tthe(indx);
                clear indx
                
                % SO
                indx=find(Gthe<=clu);
                Gso=Gso(indx);
                Tso=Tso(indx);
                clear indx
                
            case 'AnatGp2'
                
                clu=find(Map(:,2) > 5, 1,'first');
                Map=Map(clu:end, :);
                
                % Theta
                indx=find(Gthe>=clu);
                Gthe=Gthe(indx);
                Tthe=Tthe(indx);
                clear indx
                
                % SO
                indx=find(Gthe>=clu);
                Gso=Gso(indx);
                Tso=Tso(indx);
                clear indx   
                
        end
        
        BinSize = (Par.SampleRate/1000); %1 ms
        
        
        acgThe=zeros((2*HalfBins+1), size(Map,1));
        acgSo=zeros((2*HalfBins+1), size(Map,1));
        
        for n=1:size(Map,1)
            t=Tthe(Gthe==Map(n,1)); % theta
            [acg, tbin] = CCG(t, ones(size(t)), BinSize, HalfBins, Par.SampleRate,unique(ones(size(t))),'hz');
            if ~isempty(acg); acgThe(:,n)=acg; end
            
            t=Tso(Gso==Map(n,1)); % so
            [acg, tbin] = CCG(t, ones(size(t)), BinSize, HalfBins, Par.SampleRate,unique(ones(size(t))),'hz');
            if ~isempty(acg); acgSo(:,n)=acg; end
        end
        
        
        save([FileBase '.' mfilename '.' SubSet '.' HalfBins '.mat'],  'acgThe','acgSo','tbin', 'Map');
        
        
        
        
    case 'display'
        
        if FileExists([FileBase '.' mfilename '.' SubSet '.' HalfBins '.mat'])
            load([FileBase '.' mfilename '.' SubSet '.' HalfBins '.mat']);
        else
            disp('Compute first !!')
        end
        
        figure('name',['ACGs States - ' SubSet],'NumberTitle','off')
        
        for n=1:size(Map,1)
            subplot 121; bar(tbin,acgThe(:,n),'k') ; axis tight;  title(['ACG  the  ' num2str(Map(n,1))]);  ylabel('Hz')
            subplot 122; bar(tbin,acgSo(:,n),'k') ; axis tight; title(['ACG  so  ' num2str(Map(n,1))]);  ylabel('Hz')
            waitforbuttonpress
        end
        
        
end