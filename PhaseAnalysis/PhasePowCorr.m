function out = PhasePowCorr(FileBase,varargin)
%function PhasePowCorr(FileBase,fMode, ElLoc, FreqRange,State,FrThreshold,Window)
% test of Womelsdorf Science paper idea - amplitude correlations are
% highest at the preferred phase shift
[fMode, ElLoc, FreqRange,States,FrThreshold,Window] = DefaultArgs(varargin, ...
    {'compute','c',[30 150], {'REM','RUN','SWS'},10,0.25});

Par = LoadPar([FileBase '.xml']);
MinPeriod = 5; %seconds for theta periods selection
if strcmp(ElLoc,'all')
    El = [1:Par.nElecGps];
elseif isstr(ElLoc(1))
    El = find(strcmp(Par.ElecLoc,ElLoc));
else
    El = ElLoc(1);
end
if ischar(States); States = {States}; end;
nSt = length(States);

[RepCh IfTetrode Info] = RepresentChan(FileBase);
if Info.nChInShank==8 IfTetrode =1; end
nRepCh = length(RepCh);
eFs = 1250;
Window = Window*eFs;
switch fMode
    case 'compute'
        [AllRes,AllClu,Map] = LoadCluRes(FileBase,El);
        if isempty(AllRes)
            fprintf('no units in %s in this file \n',ElLoc);
            return;
        end
        AllRes = round(AllRes/16)+1;
        if ~isstr(ElLoc)
            myi = find(Map(:,2)==ElLoc(1) & Map(:,3)==ElLoc(2));
            AllRes  = AllRes(AllClu==myi);
            AllClu = ones(length(AllRes),1);
        end
        sc=1;
        %loop through States
        for s=1:nSt
            %Period = SelectStates(FileBase,States{s},MinPeriod);
            if ~FileExists([FileBase '.sts.' States{s}])
                fprintf('No state %s\n',States{s});
                continue;
            end

            Period = load([FileBase '.sts.' States{s}]);
            PerLen = round(sum(diff(Period,1,2)/1250));
            fprintf('Processing state %s of duration %d seconds\n',States{s},PerLen);

            [Res ind]=SelPerDiscr(AllRes,Period,1,1);
            Clu=AllClu(ind);

           
            fr = FiringRate(Res,Clu,[],1250);
            gfr = find(fr>FrThreshold);

            coh=[]; ph=[];powf=[];powu=[];
            % h = waitbar(0,'Please wait...');

            for ii=9%1:nRepCh
                eeg = LoadBinary([FileBase '.eeg'],RepCh(ii),Par.nChannels,Period);
                % waitbar(ii/nRepCh,h);
                for jj=1:length(gfr);
                    myRes = Res(Clu==gfr(jj));
                    if length(myRes)==0 continue; end
                    NW=3;
                    
                    [per f dfWeight] = mtptperiodogram(eeg,myRes,ones(length(myRes),1),...
                        2^(nextpow2(Window)),1250,round(Window),round(Window/2),NW,'linear',[],FreqRange,[],[],[],2*NW-1);
                   
                    [out(sc).rcorr(:,:,ii,jj) out(sc).AvR(:,:,:,ii,jj) out(sc).AvPh(:,:,:,ii,jj)] ...
                        = PhPowCorr(per,dfWeight);
                    
                end
            end
%            keyboard
            %close(h);
            out(sc).ElClu = Map(gfr,:);
            out(sc).FirRate = fr(gfr);
            out(sc).State = States{s};
            sc= sc+1;
        end %loop through states
        if nargout<1
            save([FileBase '.' mfilename '.mat'],'out');
        end

end


function [rcorr avR avPh] = PhPowCorr(per, dfWeight)
% per = nt x nf x ntapers x 2
per = per(all(dfWeight>0,2),:,:,:);
[nt,nf,ntap,nch] = size(per);
csd = sq(mean(per(:,:,:,1).*conj(per(:,:,:,2)),3));
pow=[];
pow(:,:,1) = sq(abs(mean(per(:,:,:,1).^2,3)));
pow(:,:,2) = sq(abs(mean(per(:,:,:,2).^2,3)));

%coh = csd./sqrt(pow(:,:,1).*pow(:,:,2));

ph = angle(csd);
mph = angle(mean(exp(sqrt(-1)*ph)));
mR = abs(mean(exp(sqrt(-1)*ph)));

cph = mod(dotdot(ph,'-',mph)+pi,2*pi)-pi;

%[hst cphi] = histcI(cph,linspace(-pi,pi,7));
%[hst cphi] = histcI(1-cos(cph),linspace(0,1,7));
bincent = [-pi -2*pi/3 -pi/3 0 pi/3 2*pi/3 pi];
nbins = length(bincent);
rcorr = [];
for l=1:nf
    for k=1:nbins
        dst = 1-cos(cph(:,l)-bincent(k));
        [sdst si ] = sort(dst,'ascend');
        myi = si(1:round(nt/4));
       % myamp = sq(pow(cphi(:,l)==k,l,:));
        myamp = sq(pow(myi,l,:));
        rcorr(l,k)= RankCorrelation(myamp(:,1),myamp(:,2)); 
       % scatter(myamp(:,1),myamp(:,2),20,ones(size(myamp,1),1)*k,'filled');
       subplotfit(k,nbins); plot(myamp(:,1),myamp(:,2),'.');axis tight
    end
    myamp = sq(pow(:,l,:));
    myamp(:,1)=MakeUniformDistr(myamp(:,1));
    myamp(:,2)=MakeUniformDistr(myamp(:,2));
    
    mycph = sq(cph(:,l));
    [Av Std Bins OcMap(:,:,l)] = MakeAvF(myamp,exp(sqrt(-1)*mycph),10);
    avR(:,:,l) = abs(Av);
    avPh(:,:,l) = angle(Av);
    
end

keyboard





