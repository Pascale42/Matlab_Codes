%function out = CCGSignif(Res,Clu,BinSize,HalfBins,SampleRate,Normalization,PairsToTest,Randomiz,Width,Periods)
% INPUT:
%      all CCG parameters are usual :
%      Res,Clu,BinSize,HalfBins,SampleRate,Normalization , see CCG help
%      PairsToTest - specifies the cluster numbers of pairs (per and post)
%      to test (2 column matrix)
%      Randomiz - structure that contains infromation about randomization
%      procedure with fields Type, nRand, and other params needed for the
%      procedure. Fields are:
%        Type ='shuffle' (to shuffle isi's of spike train) or 'jitter'
%        nRand - number of randomizations
%        Alpha - percentiles for confidance band (default [5 95])
%        Tau - widnow for jitter (in samples)
%      e.g Randomiz = struct('Type','shuffle','nRand',100,'Alpha',0.05)
%      Width - if not 0 will do smooth (kernel) estimates of the CCGs where
%      Width is the kernel length in samples
%      Periods - 2-column matrix of Periods to restrict all analysis to
% % OUTPUT:
%      structure with fields:
%      CCG - actual nonrandomized CCG. this and all other matrices are 
%      of dimension: #time_bins x #Pairs
%      smCCG - it's smoothed version. same smoothing as randomization
%      tbin - time bin, in msecs as usual
%      AvShufCCG - mean of CCG across randomizations
%      StdShufCCG - std of CCG across randomizations
%      ZScore = (AvShufCCG -smCCG)./StdShufCCG
%      PvalShufCCG - empirical P-value 
%      ConfBand - confidence interval for give Randomiz.Alpha 
%      GlobalPvalMax/Min - global bands p-value for minmum/maximum

function out = CCGSignif(Res,Clu,BinSize,HalfBins,SampleRate,Normalization,PairsToTest,varargin)
defRandomiz = struct('Type','jitter','nRand',100,'Tau',4,'Alpha',[5 95]);

[Randomiz, Width, Periods] = DefaultArgs(varargin,{defRandomiz,0,[]});

nPeriods = size(Periods,1);
nPairs = size(PairsToTest,1);

out.AvShufCCG = zeros(HalfBins*2+1,nPairs);
out.StdShufCCG = zeros(HalfBins*2+1,nPairs);
out.PvalShufCCG = zeros(HalfBins*2+1,nPairs);
out.CCG = zeros(HalfBins*2+1,nPairs);
out.smCCG = zeros(HalfBins*2+1,nPairs);
out.GlobalPvalMin = zeros(HalfBins*2+1,nPairs);
out.GlobalPvalMax = zeros(HalfBins*2+1,nPairs);
%out.ConfBands = zeros(HalfBins*2+1,length(Randomiz.Alpha),nPairs);
out.GlobalBandsMin = zeros(HalfBins*2+1,length(Randomiz.Alpha),nPairs);
out.GlobalBandsMax = zeros(HalfBins*2+1,length(Randomiz.Alpha),nPairs);
for p=1:nPairs
    Res1 = Res(Clu==PairsToTest(p,1));
    Res2 = Res(Clu==PairsToTest(p,2));
    if isempty(Res1) | isempty(Res2)
        warning('CCGSignif:nospks','empty spike trains');
        continue
    end
    nRes2 = length(Res2);
    if Width==0
        [T,G] = CatTrains({Res1,Res2},{1,2},Periods,1);
        [ccg out.tbin] = CCG(T,G,BinSize,HalfBins,SampleRate,[1 2],Normalization);
        smccg = sq(ccg(:,1,2));
        out.smCCG(:,p) = smccg(:);
    else
        warning('Periods are not taken into account here');
        [smccg ccg out.tbin Width] = sCCG(Res1,Res2,BinSize,HalfBins,SampleRate,Normalization,Width);
        out.smCCG(:,p) = smccg(:);
    end
    out.CCG(:,p) = sq(ccg(:,1,2));
    out.tbin = out.tbin';
    %now we keep Res1 fixed and randomize Res2
    for s=1:Randomiz.nRand
        switch Randomiz.Type
            case 'shuffle'
                if isempty(Periods) %easier case
                    isi = [0;diff(Res2)];
                    sRes2 = Res2(1)+cumsum(isi(randperm(nRes2)));
                else
                    %no time to vectorize - can be done
                    sRes2 =[];
                    for i=1:nPeriods
                        myres = Res2(Res2>Periods(i,1) & Res2<Periods(i,2));
                        if ~isempty(myres)
                            isi = [0;diff(myres)];
                            sRes2 = [sRes2; myres(1)+cumsum(isi(randperm(length(isi))))];
                        end
                    end
                end
                sRes1 = Res1;
            case 'jitter'
                ResShift = round(2*(rand(size(Res2))-0.5)*Randomiz.Tau);
                sRes2 = Res2+ResShift;
               % ResShift = round(2*(rand(size(Res1))-0.5)*Randomiz.Tau);
               % sRes1 = Res1+ResShift;
                sRes1 = Res1;
        end
        if Width==0
            [sT,sG] = CatTrains({sRes1,sRes2},{1,2},Periods,1);
            smccgr = CCG(sT,sG,BinSize,HalfBins,SampleRate,[1 2],Normalization);
            smccgr = sq(smccgr(:,1,2));
            accsmccgr(:,s)= smccgr;
        else
            [smccgr sccg tbin] = sCCG(sRes1,sRes2,BinSize,HalfBins,SampleRate,Normalization,Width);
            accsmccgr(:,s)= smccgr(:);
        end

%         out.AvShufCCG(:,p) = out.AvShufCCG(:,p) + sq(smccgr(:));
%         out.StdShufCCG(:,p) = out.StdShufCCG(:,p) + sq(smccgr(:)).^2;
%         out.PvalShufCCG(:,p) = out.PvalShufCCG(:,p) + (smccgr(:)>smccg(:));
%         
    end % end of shuffling loop
    out.AvShufCCG(:,p) = mean(accsmccgr,2);
    out.StdShufCCG(:,p) = std(accsmccgr,0,2);
    out.PvalShufCCG(:,p) = sum(accsmccgr > repmat(smccg(:),1,Randomiz.nRand),2);
    out.PointBands(:,:,p) = prctile(accsmccgr,Randomiz.Alpha,2);
    
    CCGMax = max(accsmccgr);
    CCGMin = min(accsmccgr);
   
    out.GlobalPvalMax(:,p) = sum(repmat(out.smCCG(:,p),1,Randomiz.nRand) < repmat(CCGMax(:)',2*HalfBins+1,1),2)/Randomiz.nRand;
    out.GlobalPvalMin(:,p) = sum(repmat(out.smCCG(:,p),1,Randomiz.nRand) > repmat(CCGMin(:)',2*HalfBins+1,1),2)/Randomiz.nRand;
    maxalpha = find(Randomiz.Alpha>50);
    minalpha = find(Randomiz.Alpha<50);
    if ~isempty(maxalpha)
    out.GlobalBandsMax(:,[1:length(maxalpha)],p) =  repmat(prctile(CCGMax,Randomiz.Alpha(maxalpha)),length(out.tbin),1);
    end
    if ~isempty(minalpha)
    out.GlobalBandsMin(:,[1:length(minalpha)],p) =  repmat(prctile(CCGMin,Randomiz.Alpha(minalpha)),length(out.tbin),1);
    end
end %end of pairs loop

out.PvalShufCCG = out.PvalShufCCG / Randomiz.nRand;
out.ZScore = (out.smCCG - out.AvShufCCG) ./ out.StdShufCCG;
out.GlobalBandsMax = sq(out.GlobalBandsMax);
out.GlobalBandsMin = sq(out.GlobalBandsMin);
%keyboard

