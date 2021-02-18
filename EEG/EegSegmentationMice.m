% function EegSegmentation(FileBase, Control, Overwrite, Channels, nPCA, FreqRanges,  Window, nChannels)
% performes the segmentation of the eeg by states that correspond to 
%certain power spectral features combimation
function EegSegmentation(FileBase, varargin)

esr = 1250;
[Control, Overwrite, Channels, nPCA, FreqRanges,        Window, nChannels] = DefaultArgs(varargin, ...
    {0, 1, [],           5,          [1 4; 6 9; 10 17; 35 75],  1,         []  });
saveFn = [FileBase '.eegseg.mat'];
JustDone=0;
if ~Control | ~FileExists(saveFn) % compute stuff and save
    
    if ~FileExists([FileBase '.par']) & isempty(nChannels)
        error('par file does not exist and nChannels not set');
    else
        Par = LoadPar([FileBase '.par']);
        nChannels = Par.nChannels;
    end
%     if isempty(Channels)
%         if FileExists([FileBase '.eegseg.par']);
%             Channels = load([FileBase '.eegseg.par']);
%             Channels = Channels+1;
%         elseif Control
%             system('rgd&');
%             Channels = input('So what are the channels, ah? :))\n');
%         else
%             error('Couldnt figure out the CHannels to use!!! \n');
%             return
%         end
%     end
    nChanSel =2;% length(Channels);
    SampleRate = 1e6/Par.SampleTime/16;
    Window = 2^round(log2(Window*SampleRate));
    %Eeg = readmulti([FileBase '.eeg'], nChannels, Channels);
    Eeg  = cxhipeeg(FileBase);
    Eeg = Eeg(:,[2 1]);
    % compute power in segments and get stats on it
    
    weeg = WhitenSignal(Eeg);
    [y,f,t]=mtcsglong(weeg,2^11,SampleRate,Window,[],1.5,'linear',4,[1 100]);
    %[y,f,t]=mtcsglong(eeg,2^11,1250,2^10,2^9,1.5,'linear',4,[1 100]);
    GoodF = find(f<57 | f>63); 
    
    % now do the svd on non-60Hz freq. on each channel
    for iCh = 1:nChanSel
        [u,s,v] = svd(log(y(:,GoodF,iCh)),0);
        Fet(:,:,iCh) =u(:,1:nPCA);
        sv=diag(s); sv = sv./sum(sv);
        EigVal(:,iCh) = sv;
        EigVec(:, :, iCh) = v(:,1:nPCA);
    end
    
    
    % norm Fet
    AllFet =reshape(Fet,size(Fet,1),[]);
    AllFet = round(AllFet*2^15/max(abs(AllFet(:))));
    %    AllFet = [AllFet round(t(:)*1250)];
    if 0 % do gmm clustering
        nStates = 6;
        nDim = 10;
        mix = gmm(nDim,nStates, 'full');
        opt = zeros(15,1); opt(14)=100;
        mix = gmminit(mix, AllFet,opt);
        clu = gmmclu(mix, AllFet);
    elseif 1 % or KlustaKwik
        clu = KlustaKwik(AllFet , 5);
    end
    
    
    %save stuff
    if Overwrite | ~FileExists(saveFn) 
        save(saveFn, 'y','f','t','Fet','clu','EigVal','EigVec');
    end
    JustDone=1;
end 
%%%%%%%%%%%%%%%%%%END OF COMPUTATION    
if Control % load stuff and plot
    if ~JustDone & FileExists(saveFn)  
        load(saveFn);
    end   
    
    %     % clean clu
    %     for i=4:length(clu)
    %         if clu(i)~=clu(i-1) & length(unique(clu([- 2, -1 1 2]+i)))==1
    %             clu(i)=clu(i-1);
    %         end
    %     end
    %     % plot eigenstuff
%     figure(171)
%     subplotfit(1,nPCA+1);
%     plot(log(EigVal),'.');
%     title('Eigenvalues');
%     GoodF = find(f<57 | f>63); 
%     for i=1:nPCA
%         subplotfit(i+1,nPCA+1);
%         plot(f(GoodF),abs(sq(EigVec(:,i,:))));
%         title([num2str(i) ' PCA']);
%     end
    
    % % now compute the power in particular freq. ranges.
    % for i=1:nChannels
    %     powst(:,:,:,i) = PowerPeakStats(log(yw(:,:,i)),f, FreqRanges);
    % end
    % figure
    % HistMatrix(sq(permute(powst(:,:,5,:),[2 1 4 3])),100);
    t =  t+ t(2)-t(1);
    Kluster([FileBase '.seg'],{Fet, permute(log10(y(:,1:32,:)),[3,2,1]), t ,clu},[],1250);
    
    answ = input('Do you want to save slutering changes?\n (1/0)');
    if answ
        clu = load([FileBase '.seg.clu.1']);
        nclu = clu(1); clu = clu(2:end);
    end
    uclu = unique(clu);
    nclu = max(clu);
    hclu = hist(clu,uclu);
    if length(hclu)<nclu
        hclu = [0; hclu(:)];
    end
    figure(172)
    subplotfit(1,nclu);
    bar([1:nclu],hclu)
    title('proportion of time for each state');
    
    for i=2:nclu
        subplotfit(i,nclu)
        my = find(clu==i);
        mypow = sq(mean(20*log(y(my,:,:)),1));
        mypowerr = sq(std(20*log(y(my,:,:)),0,1))/hclu(i);
        plot(f,mypow);
        hold on
        plot(f,mypow - mypowerr,'--');
        plot(f,mypow + mypowerr,'--');
        title(['cluster # ' num2str(i)]);
    end
    
    figure(173)
    CCG(t,clu,1,100,1,uclu,'count');
    ForAllSubplots('axis off');
    
    figure(174)
    AllFet =reshape(Fet,size(Fet,1),[]);
    KlustaSistData(AllFet, clu);
    
    
    
    msave('t',round(t*1250),'w');
    msave('c',[nclu; clu(:)],'w');
    
    system('rgd&');
    %     %show clusters 
    %     Kluster([FileBase '.seg'],{Fet, permute(log10(y(:,1:32,:)),[3,2,1]), t ,clu},[],1250);
    
    
    %now join times in segments
    % add zeros at the ends to fix the borders
    Change = diff([ 0; clu(:); 0]);
    ChangeInd  = find(abs(Change)>0.5);
    ChangeISI = diff(ChangeInd);
    Seg = [ChangeInd(1:end-1) ChangeInd(1:end-1)+ChangeISI(:)-1];
    SegClu = clu(ChangeInd(1:end-1));
    SegLen = ChangeISI(:);
    hSegClu = hist(SegClu,uclu);
    if length(hSegClu)<nclu
        hSegClu = [0; hSegClu(:)];
    end
    
    Labels = {'SWS', 'IMM', 'REM', 'RUN' , 'AWK', 'HVS' , 'ART' ,'DES'};
    % make browsing sequence for each cluster
    for i=2:nclu
        myclu = find(SegClu==i);
        myt = t(round(mean(Seg(myclu,:),2)));
        [dummy ind] = sort(SegLen(myclu));
        nb = min(20, hSegClu(i));
        % take half from the top
        settop = ind(end-round(nb/2):end);
        % and another half randomly from the bottom half
        setbottom = round( rand(round(nb/2),1) * hSegClu(i)/2)+1;
        set = [settop(:) ; setbottom(:)];
        b = round(myt(set(:))*1250);
        msave(['b' num2str(i)],b(:),'w');
        msave('b' ,b(:),'w');
        fprintf('assign label to the cluster %d\n', i);
        fprintf('chose:  SWS (1), IMM (2), REM (3), RUN (4), AWK(5), HVS (6), ART (7), DES (8) \n');
        CluLabel(i) = input(' ----> ');
    end
    CluLabel(1)=7;
    fprintf(' you assigned: \n');
    for i=1:nclu
        fprintf('Cluster %d ---> %s\n', i, Labels{CluLabel(i)});
    end
    answ = 1;
    while answ
        answ = input('Do you want to reasign ? (1/0)')
        if answ 
            whichclu = input('what cluster? \n');
            fprintf('assign label to the cluster %d\n', whichclu);
            fprintf('choose:  SWS (1), IMM (2), REM (3), RUN (4), HVS (5), ART (6), DES (7) \n');
            CluLabel(i) = input(' ----> ');
        end
    end
    CluLabelStr = Labels(CluLabel);
    SegClu = CluLabel(SegClu(:));
    
    
    if 0    
        %Set the running borders
        Par = LoadPar([FileBase '.par']);
        SyncChannel = Par.nChannels; % assume for now;
        Sync = readsinglech([FileBase '.eeg'],Par.nChannels,SyncChannel);
    end
    %     figure(175)
    %     subplot(3,1, [1 2];
    %     stairs(Seg(:,1)*(t(2)-t(1)),CluLabel(SegClu));
    %     ylim([0 8]);
    % %     ylables =get(gca, ;
    % %     for i=1:6
    % %         ylabels(i,:) = Labels{i};
    % %     end
    % %     ylabels = ['  ';ylabels; '   '];
    % %     set(gca,'YTickLabel', ylabels)
    %     subplot(3,1,3);
    %     plot([1:length(Sync)], Sync);
    %     [runbound,dummy] = ginput(2);
    
    if 0 % skip for now.
        % merging rules
        % SWS - HVS, SWS - IMM, ART-AWK,/RUN, DES-AWK/RUN/REM 
        nSegClu = length(Labels);
        %rename the segments
        SegClu = SegClu(:);
        Merge1= [ 1 6; 1 2; 7 4; 7 5;  8 4; 8 5; 8 3]; 
        Merge1 = [ Merge1 ; fliplr(Merge1)];
        
        
        figure(175)
        stairs(Seg(:,1)*(t(2)-t(1)),SegClu);
        ylim([0 9]);
        fprintf('Now will merge the single eopchs (of length 1) (hit key) \n');
        % now sort out the segments according to labels
        %length one segments
        
        WasMerge = 2;
        while WasMerge>0
            
            for i=[-1 1]
                Seg1Ind = find(SegLen<3 );
                Seg1NeighbInd = Seg1Ind + i;
                NeighbLen  = SegLen(Seg1NeighbInd);
                NeighbClu = SegClu(Seg1NeighbInd);
                % wjat to merge
                IfMerge = ismember([NeighbClu SegClu(Seg1Ind)], Merge1,'rows') & (NeighbLen>2);
                if sum(IfMerge)==0 
                    WasMerge = WasMerge -1 ;     
                    break;
                end
                MergeInd = Seg1Ind(find(IfMerge)); % in Seg indexing
                % now we merge
                Remaining = setdiff([1:length(SegLen)], MergeInd);
                SegClu = SegClu(Remaining);
                if i==-1
                    Seg(MergeInd-1, 2) = Seg(MergeInd-1, 2)+1;
                    SegLen(MergeInd-1) = SegLen(MergeInd-1)+1;
                else
                    Seg(MergeInd+1, 1) = Seg(MergeInd-1, 2) -1;
                    SegLen(MergeInd+1) = SegLen(MergeInd+1)+1;
                end
                SegLen = SegLen(Remaining);
                Seg = Seg(Remaining,:);
            end
            hold on
            stairs(Seg(:,1)*(t(2)-t(1)),SegClu, 'r');
            hold off
            pause
            stairs(Seg(:,1)*(t(2)-t(1)),SegClu, 'b');
        end
        
    end % end of merge section
    close all
    plot(Seg(:,1),SegClu,'.','MarkerSize',3)

    Seg = round(Seg*(t(2)-t(1))*1250);
    msave([FileBase '.states.res'], Seg, 'w');
    msave([FileBase '.states.len'], SegLen(:),'w');
    msave([FileBase '.states.clu'], SegClu(:),'w');
end



return

