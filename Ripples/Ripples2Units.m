%
% function OutArgs = Ripples2Units(FileBase,fMode, varargin :: Normalization, SubSet, SpkGps, CluLocOnly)
%
% fMode :: 'compute', 'display'
% Normalization :: 1 or 0; normalizes correlation on ripples duration
% SubSet : Name it as you wish or 'All' ;
% SpkGps : enter all the SpikeGroup numbers, ex for AnatGp2: [6:10] or [6 7 8 9 10] (default = [])
% CluLocOnly : run the function only on the clusters recorded by this list of channels (from 1) (default = [])
% HalfWin : Half duration of the correlation window, in ms
% jscale : size of jitter window (in ms) for statistics (when Normalization = -1)
%
% plots show neuronal firing  +/- 200ms around ripple time
% If Normalization, ripple duration is [-1 1], otherwise the mean duration is shown by a dashed orange line
% If Normalization = -1; Computes statistics on each ripple/neuron pair


% Pascale heavily modified : 30 avril 2015, 5 mai 2015 and 29 juillet 2015


function OutArgs = Ripples2Units(FileBase,fMode, varargin)
[Normalization, SubSet, SpkGps, CluLocOnly, HalfWin, jscale] =  ...
    DefaultArgs(varargin,{0,'All',[],[],200,5});

switch fMode
    case 'compute'
        Par = LoadPar([FileBase '.xml']);
        
        if ~FileExists([FileBase '.spw.mat']) && ~FileExists([FileBase '.DetectRipples.mat'])
            fprintf('no Ripples detected');
            return
        end
        
        %         if FileExists([FileBase '.spw.mat'])
        %              load([FileBase '.spw.mat']);
        %         end
        %
        %        if   FileExists([FileBase '.DetectRipples.mat'])
        %            load([FileBase '.DetectRipples.mat'])
        %            Rips=OutArgs.Ripples;
        %        end
        
        if FileExists([FileBase '.SWS.HPC_RIPPLES_STE.mat'])
            Rips=load([FileBase '.SWS.HPC_RIPPLES_STE.mat']);
            Rips=Rips.HPC.st_HFOInfo.m_EvtLims;
            SWS=load([FileBase '.sts.SWS']);
            Rips=Rips+SWS(1,1);
            Rips(365:end,:)=Rips(365:end,:)+3144776;
            
        else
            return
        end
        
        
        %%% Loading units
        
        if strcmp(SubSet, 'All');
            if ~FileExists([FileBase '.CluRes.mat'])
                [T,G,Map]=LoadCluRes(FileBase);
            else
                load([FileBase '.CluRes.mat']);
            end
        else
            [T,G,Map]=LoadCluRes(FileBase, SpkGps,0, CluLocOnly);
        end
        
        T =  round(T*Par.lfpSampleRate/Par.SampleRate);
        
        
        %%% Ripples 2 Units
        
        if Normalization
            uclu=unique(G);
            OutArgs.ccg=[];
            
            for clu=1:length(uclu) % for each neuron
                t=T(G==uclu(clu));
                [tsOffsets, ts1idx, ts2idx] = crosscorrelogram(Rips.t/Par.lfpSampleRate, t/Par.lfpSampleRate, [-0.2 0.2]);
                for n=1:length(tsOffsets)
                    tsOffsetsN(n)= tsOffsets(n)/Rips.len(ts1idx(n));
                end
                [ccg, OutArgs.tbin]=hist(tsOffsetsN, 401);
                OutArgs.ccg=[OutArgs.ccg ; ccg];
            end
            
            OutArgs.ccg=OutArgs.ccg';
            OutArgs.norm=1;
            
        end
        
        if Normalization == 0
            [ccg, OutArgs.tbin] = Trains2CCG({Rips(:,1),T},{1,G},1,125,Par.lfpSampleRate,'count');
            %  ccg = squeeze(ccg(:,1,2:end));
            OutArgs.ccg = squeeze(ccg(:,1,2:end));
            OutArgs.norm=0;
        end
        
        if Normalization == -1
            %            [ccg, OutArgs.tbin] = Trains2CCG({Rips.t,T},{1,G},1,125,Par.lfpSampleRate,'scale');
            %           %  ccg = squeeze(ccg(:,1,2:end));
            %             OutArgs.ccg = squeeze(ccg(:,1,2:end));
            OutArgs.norm=-1;
            
            %%% Jitter
            njitter = 1000;
            jscale=50;
            
            
            for n=Map(:,1)'
                
                [ccg, tbin] = CCG([Rips.t;T(G==n)],[ones(size(Rips.t));2*ones(size(T(G==n)))], 1, ...
                    (200*Par.lfpSampleRate/1000), Par.lfpSampleRate,[1,2],'count');
                OutArgs.Ccg(:,n)=squeeze(ccg(:,1,2));
                OutArgs.tbin=tbin;
                disp(['compute neuron ' num2str(n)])
                
                
                %  CCG for jittering data
                for x=1:njitter
                    T_jitter = T(G==n) + round(2*round(Par.lfpSampleRate/1000*jscale)*rand(size(T(G==n)))-1 ...
                        *round(Par.lfpSampleRate/1000)*jscale);
                    [jccg, tJ] = CCG([Rips.t;T_jitter],[ones(size(Rips.t));2*ones(size(T_jitter))], 1, ...
                        (200*Par.lfpSampleRate/1000), Par.lfpSampleRate,[1,2],'count');
                    ccgj(:,x)=jccg(:,1,2);
                end
                %  Computes the pointwise line
                ccgjm(:,n)  = mean(ccgj,2);
                ccgjmax(:,n)=ccgjm(n) + 2*std(ccgj,[],2);
                ccgjmin(:,n)=ccgjm(n) - 2*std(ccgj,[],2);
            end
            
            OutArgs.J.ccgjm=ccgjm;
            OutArgs.J.ccgjptMax=ccgjmax;
            OutArgs.J.ccgjptMin=ccgjmin;
            OutArgs.J.njitter=njitter;
            OutArgs.J.jscale=jscale;
            
        end
        
        OutArgs.map = Map;
        OutArgs.RipLen=mean(Rips(:,2)-Rips(:,1));
        
        %%% Std firing
        
        % TO DO
        
        
        
        save([FileBase '.' mfilename '.' SubSet '.mat'],'OutArgs');
        
    case 'display'
        load([FileBase '.' mfilename '.' SubSet '.mat'],'OutArgs');
        if ~isempty(OutArgs)
            
            figure('name',[mfilename ' - ' FileBase ' - ' SubSet],'NumberTitle','off')
            ncells =size(OutArgs.ccg,2);
            for n=1:ncells
                subplotfit(n,ncells);
                if OutArgs.norm == 0
                    bar(OutArgs.tbin, squeeze(OutArgs.ccg(:,n))', 'FaceColor',[0.4 .4 .4],'EdgeColor',[0.4 .4 .4],'LineWidth',0.5);
                    hold on
                    line([0 0], [0, max(OutArgs.ccg(:,n))], 'Color', [1 0 0], 'LineWidth', 2)
                    axis tight
                    title(num2str(OutArgs.map(n,2:3)));
                    %                     line([-(OutArgs.RipLen/2) -(OutArgs.RipLen/2)], [0 max(OutArgs.ccg(:,n))],'linestyle','--','Color',[1 0.5 0],'LineWidth',0.5)
                    %                     line([OutArgs.RipLen/2 OutArgs.RipLen/2], [0 max(OutArgs.ccg(:,n))],'linestyle','--','Color',[1 0.5 0],'LineWidth',0.5)
                end
                
                if OutArgs.norm == -1
                    bar(OutArgs.tbin, squeeze(OutArgs.Ccg(:,n))', 'FaceColor',[0.4 .4 .4],'EdgeColor',[0.4 .4 .4],'LineWidth',0.5);
                    hold on
                    line([0 0], [0, max(OutArgs.Ccg(:,n))], 'Color', [1 0 0], 'LineWidth', 2)
                    axis tight
                    title(num2str(OutArgs.map(n,2:3)));
                    % Significance
                    line(OutArgs.tbin,OutArgs.J.ccgjm,'linestyle','--','color','b','LineWidth',1.5)
                    line(OutArgs.tbin,OutArgs.J.ccgjptMax,'linestyle','--','color','r','LineWidth',1.5)
                    line(OutArgs.tbin,OutArgs.J.ccgjptMin,'linestyle','--','color','r','LineWidth',1.5)
                    set(gca,'XLim',[min(OutArgs.tbin),max(OutArgs.tbin)]); axis tight
                end
                
                if OutArgs.norm == 1
                    bar(OutArgs.tbin.*1000, squeeze(OutArgs.ccg(:,n))', 'FaceColor',[0.4 .4 .4],'EdgeColor',[0.4 .4 .4],'LineWidth',0.5);
                    hold on
                    line([0 0], [0, max(OutArgs.ccg(:,n))], 'Color', [1 0 0], 'LineWidth', 2)
                    axis tight
                    title(num2str(OutArgs.map(n,2:3)));
                    xlim([-2000 2000])
                    %                     NON!!!!!!
                    %                     line([-1 -1], [0 max(OutArgs.ccg(:,n))],'linestyle','--','Color',[1 0.5 0],'LineWidth',0.5)
                    %                     line([1 1], [0 max(OutArgs.ccg(:,n))],'linestyle','--','Color',[1 0.5 0],'LineWidth',0.5)
                    %                     xlim([-4 4])
                end
            end
            
        else
            return
        end
end
