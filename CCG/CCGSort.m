function Count = CCGSort(Pairs, Res, Clu, SortInd, SortValue)
nbins = 10;

MyPairs = Pairs(Clu(Pairs(:,1))==1 & Clu(Pairs(:,2))==2,:);
BinInd = round(diff(Res(MyPairs),1,2)/CCGBinSize) + HalfBins + 1;

SortBins = linspace(min(SortValue),max(SortValue),nbins+1);

[ValHist ValInd] = histcI(SortValue, SortBins);
gbins = find(ValInd>0);
Count= Accumulate([BinInd(gbins), ValInd(gbins)], 1, [2*HalfBins+1, nbins]);

% n1 = histcI(ThAmp(Res(Clu==1)), ampbins);
% n2 = histcI(ThAmp(Res(Clu==2)), ampbins);
% Expected = n1.*n2/(2*HalfBins+1);
%out(cnt).ccg_amp = CountAmp; %(CountAmp+Eps) ./(repmat(Expected',2*HalfBins+1,1)+Eps);


