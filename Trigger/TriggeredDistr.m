%function [HistMatr, Time, HistBins] =TriggeredDistr(Trace, nBefore, nAfter, T, G, nBins, IfNorm, AmpSize, Sigma, Fs)
%
% computes triggered distribution from Trace at the times given by T.
% Trace may be 2D, in which case columns are averaged separately.
% The output will then be of the form HistMat(Bins,Time, Column)
% nBefore and nAfter give the number of samples before and after
% to use.
%
% G is a group label for the trigger points T.  In this case the
% output will be HistMat(Bins, Time, Column, Group)
% Time and HistBins give the coordinates for matrix plot
% use ImageMatrix(Time, HistBins, HistMat) to visualize

function [HistMatr, Time , HistBins] = TriggeredDist(Trace, nBefore, nAfter, T, varargin)

[G, nBins, IfNorm, AmpSize, Sigma, Fs] = DefaultArgs(varargin, { ones(length(T), 1), 20, 1, 1.6, 0.005,1250});

nColumns = size(Trace,2);
nSamples = nBefore + nAfter + 1;
nGroups = max(G);
maxTime = size(Trace, 1);

%normalize if IfNorm
if IfNorm 
    Trace = NormMatrix(Trace);
end
if ~IfNorm 
    HistBinsAll = zeros(nBins,nColumns,nGroups);
else
    HistBinsAll = repmat((linspace(-AmpSize, AmpSize, nBins))',[1, nColumns,nGroups]);
end
BlockSize = floor(2000000/nSamples); % memory saving parameter
%XTriggered = zeros(nTriggers,nSamples);
HistMatr = zeros(nBins,nSamples,nColumns, nGroups);
BinCnt = zeros(nBins,nSamples,nColumns, nGroups);
for grp = 1:nGroups

    Sum = zeros(nSamples, nColumns);
    SumSq = zeros(nSamples, nColumns);
    MyTriggers = find(G==grp & T > nBefore & T <= maxTime-nAfter);
    nTriggers = length(MyTriggers);
    
    % go through triggers in groups of BlockSize to save memory
    for Block = 1:ceil(nTriggers/BlockSize)
        BlockTriggers = MyTriggers(1+(Block-1)*BlockSize:min(Block*BlockSize,nTriggers));
        nBlockTriggers = length(BlockTriggers);
        
        TimeOffsets = repmat(-nBefore:nAfter, nBlockTriggers, 1);
        TimeCenters = repmat(T(BlockTriggers), 1, nSamples);
        TimeIndex = TimeOffsets + TimeCenters;
        
        Waves = Trace(TimeIndex,:);
        Waves = reshape(Waves, [nBlockTriggers, nSamples, nColumns]);
        for i=1:nColumns

            if Block==1& ~IfNorm
                % find hist and bins from first Block if not normalized
                [BlockHist, HistBins] =hist(squeeze(Waves(:,:,i)) ,nBins);
                BinCnt(:,:,i,grp) = BinCnt(:,:,i,grp)+ones(nBins,nSamples)*nBlockTriggers;
                HistBinsAll(:,i,grp) = HistBins;
            else
                HistBins = squeeze(HistBinsAll(:,i,grp));
                delBin = HistBins(2,:) - HistBins(1,:);
                BinEdges = HistBins + repmat(delBin/2, nBins,1);
                BinEdges = [HistBins(1,:)-delBin/2; BinEdges];
                [BlockHist ind]= histcI(squeeze(Waves(:,:,i)) , BinEdges);
                BinCnt(:,:,i,grp) = BinCnt(:,:,i,grp) + repmat(sum((ind>0),1),nBins,1);
             end

            HistMatr(:,:,i,grp) = HistMatr(:,:,i,grp) + BlockHist;
        end
    end
    
end
HistMatr = HistMatr ./ BinCnt ;
if (nargout>1 & Sigma~=0)
     HistMatr = SmoothMatrix2(HistMatr, 0, Sigma);   
 end
Time = [-nBefore:nAfter];
HistBins = HistBinsAll;
if (nargout<1)
    figure
    HistMatr = SmoothMatrix2(HistMatr, 0, Sigma);    
    ImageMatrix(Time*1000/Fs, HistBins, HistMatr);
end

