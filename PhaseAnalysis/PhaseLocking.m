function out = PhaseLocking(FileBase,varargin)
%function PhaseLocking(FileBase,fMode, ElLoc, FreqRange,State,FrThreshold)
% computes unit phase modulation of power at various freq. on each channel for each cell.
[fMode, ElLoc, FreqRange,States,FrThreshold] = DefaultArgs(varargin, ...
    {'compute','c',[], {'REM','RUN'},3});

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
bs=2;
n=51;
d1 = log(1)/log(bs);
d2 = log(150)/log(bs);
y = (bs).^ [d1+(0:n-2)*(d2-d1)/(floor(n)-1), d2];
FreqRange(:,1) = ((y(1:end-1)+y(2:end))/2)';
FreqRange(:,2) = diff(y(:));

switch fMode
    case 'compute'
        [AllRes,AllClu,Map] = LoadCluRes(FileBase,El);
        if isempty(AllRes)
            fprintf('no units in %s in this file \n',ElLoc);
            return;
        end
        AllRes = round(AllRes/16)+1;
        if ~isstr(ElLoc)
            myi = find(elclu(:,1)==ElLoc(1) & elclu(:,2)==ElLoc(2));
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
            [uClu dummy Clu] = unique(Clu);
            out(sc).Map = Map(uClu,2:3);
            fr = FiringRate(Res,Clu,[],1250);
 %           gfr = find(fr>FrThreshold);
            nf = size(FreqRange,1);
            
            for ii=1:nRepCh
                eeg = LoadBinary([FileBase '.eeg'],RepCh(ii),Par.nChannels,Period)';
                for jj=1:size(FreqRange)
                    %WinLength = 2^nexpow2(5*1250/FreqRange(1));
                    myFreqRange = FreqRange(jj,1)+[-0.5 0.5]*FreqRange(jj,2);
                    %feeg = MTFilter(eeg,FreqRange(jj,:), 1250,WinLength
                    feeg = ButFilter(eeg,2,myFreqRange/625,'bandpass');
                    hlb = hilbert(feeg);
                    ph = angle(hlb);
                    
                    rt = RayleighTest(ph(Res),Clu,max(Clu));
                    
                    out(sc).rt(ii,jj).logZ = rt.logZ;
                    out(sc).rt(ii,jj).n = rt.n;
                    out(sc).rt(ii,jj).th0 = rt.th0;
                    out(sc).rt(ii,jj).k = rt.k;
                end
            end
            out(sc).FirRate = fr;
            out(sc).State = States{s};
            sc= sc+1;
        end %loop through states
        save([FileBase '.' mfilename '.mat'],'out');
        
end
        