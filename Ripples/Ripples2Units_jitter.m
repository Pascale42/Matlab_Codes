%
% function OutArgs = Ripples2Units_jitter(FileBase,fMode, varargin :: SubSet, SpkGps, CluLocOnly, jscale,njitter,alph )
%
% fMode :: 'compute', 'display'
% SubSet : Name it as you wish or 'All' ;
% SpkGps : enter all the SpikeGroup numbers, ex for AnatGp2: [6:10] or [6 7 8 9 10] (default = [])
% CluLocOnly : run the function only on the clusters recorded by this list of channels (from 1) (default = [])
%  jscale :: jittering scale, in 'ms' (default = 5)
%  njitter :: # of jittering  (default = 1000)
%  alph :: significant level (default = 0.01)

% Heavily modified  30 avril 2015  & 5 mai 2015


function OutArgs = Ripples2Units_jitter(FileBase,fMode, varargin)
[SubSet, SpkGps, CluLocOnly, jscale,njitter,alph] =  ...
    DefaultArgs(varargin,{'All',[],[],5,1000,0.01});

switch fMode
    case 'compute'
        Par = LoadPar([FileBase '.xml']);

        if ~FileExists([FileBase '.spw.mat'])
            fprintf('no spws');
            return;
        end
        load([FileBase '.spw.mat']);
        rips = Rips.t;

       %%% Loading units
        
        if strcmp(SubSet, 'All');
            if ~FileExists([FileBase '.CluRes.mat'])
                [T,G,Map, Par]=LoadCluRes(FileBase);
            else
                load([FileBase '.CluRes.mat']);
            end
        else
            [T,G,Map,Par]=LoadCluRes(FileBase, SpkGps,0, CluLocOnly);
        end
        
        
        nClu = max(G);
        T =  max(1,round(T*Par.lfpSampleRate/Par.SampleRate)); 

           %%% Ripples 2 Units
           
           Bins=1;
           HalfBins=200;
        
        [ccg, OutArgs.tbin] = Trains2CCG({rips,T},{1,G},Bins,HalfBins,Par.lfpSampleRate,'scale');
        OutArgs.ccg = squeeze(ccg(:,1,2:end)); 
        OutArgs.map = Map;
        clear ccg Map
        
        %%% Jitter and significance
        

        %  Trains2CCG for jittering data
        
        for m=1:nClu
            
            indx= G == m;
            g=G(indx);
            t=T(indx);
            
            for n=1:njitter
                T_jitter = t + 2*(20*jscale)*rand(size(t))-1*20*jscale;
                [ccg, tJ] = Trains2CCG({rips,T_jitter},{1,g},Bins,HalfBins,Par.lfpSampleRate,'scale');
                ccgj(:,n)=ccg(:,1,2);
                ccgjmax(n)=max(ccgj(:,n));
                ccgjmin(n)=min(ccgj(:,n));
            end
            
            %  Computes the pointwise line
            signifpoint = njitter*alph;
            for i=1:length(tJ)
                sortjitterDescend  = sort(ccgj(i,:),'descend');
                sortjitterAscend   = sort(ccgj(i,:),'ascend');
                ccgjptMax(i) = sortjitterDescend(signifpoint);
                ccgjptMin(i) = sortjitterAscend(signifpoint);
            end
            
            %         %  Compute the global line
            %         sortgbDescend   = sort(ccgjmax,'descend');
            %         sortgbAscend    = sort(ccgjmin,'ascend');
            %         ccgjgbMax  = sortgbDescend(signifpoint)*ones(size(tJ));
            %         ccgjgbMin  = sortgbAscend(signifpoint)*ones(size(tJ));
            %
            ccgjm  = mean(ccgj,2);
            
            OutArgs.ccgJ.m(:,m)=ccgjm;
            OutArgs.ccgJ.ptMax(m,:)=ccgjptMax;
            %         OutArgs.ccgJ.gbMax=ccgjgbMax;
            OutArgs. ccgJ.ptMin(m,:)=ccgjptMin;
            %         OutArgs.ccgJ.gbMin=ccgjgbMin;
            
        end

        save([FileBase '.' mfilename '.' SubSet '.mat'], 'OutArgs', 'alph', 'njitter', 'jscale');

    case 'display'
        load([FileBase '.' mfilename '.' SubSet '.mat']);
        if ~isempty(OutArgs)

             figure('name',[mfilename ' - ' SubSet],'NumberTitle','off')
            ncells =size(OutArgs.ccg,2);
            for n=1:ncells
                subplotfit(n,ncells);
                bar(OutArgs.tbin, squeeze(OutArgs.ccg(:,n))', 'FaceColor',[0.4 .4 .4],'EdgeColor',[0.4 .4 .4],'LineWidth',0.5 ); hold on
                line([0 0], [0, max(OutArgs.ccg(:,n))], 'Color', [1 0 0], 'LineWidth', 2)
                line(OutArgs.tbin, repmat(mean(OutArgs.ccgJ.m(:,n)), length(OutArgs.ccgJ.m(:,n)),1) ,'linestyle','--','Color',[1 0.5 0],'LineWidth',0.5)
                line(OutArgs.tbin,repmat(min(OutArgs.ccgJ.ptMax(n,:)), length(OutArgs.ccgJ.ptMax(n,:)),1),'linestyle','--','Color',[0.5 1 0],'LineWidth',1)
                line(OutArgs.tbin,repmat(max(OutArgs.ccgJ.ptMin(n,:)), length(OutArgs.ccgJ.ptMin(n,:)),1),'linestyle','--','Color',[0.5 1 0],'LineWidth',1)
                
                
%                 line(OutArgs.tbin, OutArgs.ccgJ.m(:,n),'linestyle','-','Color',[1 0.5 0],'LineWidth',0.5)
%                 line(OutArgs.tbin,OutArgs.ccgJ.ptMax(n,:),'linestyle','-','Color',[1	0.63	0.48],'LineWidth',0.5)
%                 line(OutArgs.tbin,OutArgs.ccgJ.ptMin(n,:),'linestyle','-','Color',[1	0.63	0.48],'LineWidth',0.5)
%                 line(OutArgs.tbin,OutArgs.ccgJ.gbMax,'linestyle','--','Color',[0.5 1 0],'LineWidth',0.5)
%                 line(OutArgs.tbin,OutArgs.ccgJ.gbMin,'linestyle','--','Color',[0.5 1 0],'LineWidth',0.5)
%                 
                axis tight
                title(['Neuron ' num2str(OutArgs.map(n,2:3)) ' vs Ripples']);
            end
        else
            return
        end
end
