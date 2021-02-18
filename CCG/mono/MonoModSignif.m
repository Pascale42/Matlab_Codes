%function out = MonoModSignif(ThPh,Res,Clu,MaxLag,PairsToTest,Randomiz,SampleRate)
% INPUT:
    % all CCG parameters are usual :
    % ThPh is a continuous vector of phases
    % Res,Clu
    % MaxLag - till what lag to compute the cumulative stats (1 ms
    % incriment)
    % PairsToTest - specifies the cluster numbers of pairs (per and post)
    % to test (2 column matrix)
    % Randomiz - structure that contains infromation about randomization
    % procedure with fields nRand, and other params needed for the
    % procedure
    %   nRand - number of randomizations
    %   nSubs - number of uniform subsamples in phase histogram
    %   Alpha - percentiles for confidance band (default [5 95])
    %   Tau - widnow for jitter (in samples)
    % Width - if not 0 will do smooth (kernel) estimates of the CCGs where
    % Width is the kernel length in samples
    % Periods - 2-column matrix of Periods to restrict all analysis to
% OUTPUT:
function out = MonoModSignif(ThPh,Res,Clu,MaxLag,PairsToTest,varargin)
defRandomiz = struct('nRand',100,'nSubs',20,'Tau',4,'Alpha',[5 95]);

[Randomiz,SampleRate] = DefaultArgs(varargin,{defRandomiz,20000});

nPairs = size(PairsToTest,1);

for p=1:nPairs
   
    Res1 = Res(Clu==PairsToTest(p,1));
    Res2 = Res(Clu==PairsToTest(p,2));
    if isempty(Res1) | isempty(Res2)
        warning('CCGSignif:nospks','empty spike trains');
        continue
    end
    nRes2 = length(Res2);
    try
        %actual data
        out(p).lmod = TlagMod(ThPh,Res1,Res2,MaxLag,SampleRate,Randomiz.nSubs);
        
        if Randomiz.nRand>0
            %now we keep Res1 fixed and jitter Res2
            for s=1:Randomiz.nRand
                ResShift = round(2*(rand(size(Res2))-0.5)*Randomiz.Tau);
                sRes2 = Res2+ResShift;
                lmod = TlagMod(ThPh,Res1,Res2,MaxLag,SampleRate,Randomiz.nSubs);
                out(p).lmodjit.lag = lmod.lag;
                out(p).lmodjit.r(:,s) = lmod.r(:);
                out(p).lmodjit.th(:,s) = lmod.th(:);
            end % end of shuffling loop
        end
        out(p).Randomiz = Randomiz;
    catch
        out(p).lmod = NaN;
    end
end %end of pairs loop


return
function lmod = TlagMod(ThPh,Res1,Res2,MaxLag,SampleRate,nSubs)
  
[BinInd1 Cnt1 nBins1] = BinSpk(ThPh(round(Res1/16)));
[BinInd2 Cnt2 nBins2] = BinSpk(ThPh(round(Res2/16)));
r=NaN*zeros(MaxLag*2-2,nSubs);
th0=r;n=r; 
hst=NaN*zeros(MaxLag*2-2,nSubs,8);
for k=1:nSubs
    %now subsample Res1 and Res2 such that ThPh(Res) is uniform        
    UniInd1 = UniformSubset(BinInd1,Cnt1,nBins1);
    UniInd2 = UniformSubset(BinInd2,Cnt2,nBins2);
    [T,G] = CatTrains({Res1(UniInd1),Res2(UniInd2)},{1,2});
   % [T,G] = CatTrains({Res1,Res2},{1,2});
    
    [T si] = sort(T);
    G= G(si);
    [ccg tbin pairs] = CCG(T,G,floor(SampleRate/1000),MaxLag,SampleRate,[1 2],'count');
    if isempty(pairs) continue; end
    myPairs = pairs(G(pairs(:,1))==1 & G(pairs(:,2))==2,:);%???
    tLag = diff(T(myPairs),1,2)/SampleRate*1000;
    myPh = ThPh(round(T(myPairs(:,1))/16));
    lagEdge = [1:MaxLag];
    lc=0;
    for l=length(lagEdge):-1:2
        myInd = find(tLag>-lagEdge(l) & tLag<-lagEdge(1));
        lc=lc+1;
        lmod.lag(lc) = -lagEdge(l);
        if length(myInd)<2 continue; end
        m = mean(exp(i*myPh(myInd)));
        
        th0(lc,k) = angle(m);
        r(lc,k) = abs(m);
        n(lc,k) = length(myInd);
        hst(lc,k,:) = hist(myPh(myInd),linspace(-pi,pi,8));
    end
    for l=2:length(lagEdge)
        myInd = find(tLag>lagEdge(1) & tLag<lagEdge(l));
        lc=lc+1;
        
        lmod.lag(lc) = lagEdge(l);

        if length(myInd)<2 continue; end
        m = mean(exp(i*myPh(myInd)));
        
        th0(lc,k) = angle(m);
        r(lc,k) = abs(m);
        n(lc,k) = length(myInd);
        hst(lc,k,:) = hist(myPh(myInd),linspace(-pi,pi,8));
    end
end
m = nanmean(r.*exp(i*th0),2);
lmod.r = sq(abs(m));
lmod.th = sq(angle(m));
lmod.mhst = sq(nanmean(hst,2));
    
    

return


function [BinInd Cnt nBins] = BinSpk(Ph)

nBins = 16;
while 1
    binEdge = linspace(-pi-eps,pi+eps,nBins+1);
    [Cnt BinInd] = histcI(Ph,binEdge);
    mcnt = min(Cnt);
    if mcnt>0 
        break;
    else
        nBins=nBins-1;
    end
    if nBins==2
        warning('too few spikes to bin the phase');
        BinInd = 0;
        return
    end
end
return

function UniInd = UniformSubset(BinInd,Cnt,nBins)
UniInd = [];
mCnt = min(Cnt);
for b=1:nBins 
    if Cnt(b)>mCnt
        UniInd =[UniInd; randsample( find(BinInd==b), mCnt)];
    else
        UniInd =[UniInd; find(BinInd==b)];
    end
end

return