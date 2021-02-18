function [smCCG, tbin] = SmoothCCG(Res,Clu,Pairs,Bins,smN, SampleRate)
%function smCCG = SmoothCCG(Res,Clu,Pairs,Bins,smN,SampleRate)
nClu = max(Clu);
GSubset = unique(Clu);
nSpk = Accumulate(Clu,1,nClu);
%smN = round(0.1*median(nSpk));
if length(Bins)==1
    maxt = max(diff(Res(Pairs),1,2))/SampleRate*1000;
    tbin = linspace(-maxt,maxt,Bins);
else
    tbin=Bins;
end
for i1=1:length(GSubset)
    Cell1 = GSubset(i1);
    for i2=i1:length(GSubset)
        Cell2=GSubset(i2);
        
		MyPairs = Pairs(Clu(Pairs(:,1))==Cell1 & Clu(Pairs(:,2))==Cell2,:);
        nVal = size(MyPairs,1);
        if nVal>5
            
            dT = diff(Res(MyPairs),1,2)/SampleRate*1000;
%             mysmN = round(smN*nVal);
%             if mysmN<=3 mysmN=3; end
            [f,x]= smpdf(dT, smN);
            [ux, ui]=unique(x);

            smCCG(:,i1,i2) = interp1(ux,f(ui),tbin,'cubic','extrap');
        end
    end
end

if nargout<1
%     figure
    BarMatrix(tbin, smCCG);
end
