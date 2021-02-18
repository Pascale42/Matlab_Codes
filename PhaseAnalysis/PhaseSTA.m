function out = PhaseSTA(FileBase,varargin)
%function PhaseSTA(FileBase,fMode, ElLoc, State,Window)
[fMode, ElLoc,States,Window] = DefaultArgs(varargin, ...
    {'compute','c', {'REM','RUN'},25});

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

[RepCh ifTetr Info] = RepresentChan(Par);
nRepCh = length(RepCh);
eFs = 1250;

switch fMode
    case 'compute'
        [AllRes,AllClu,Map] = LoadCluRes(FileBase,El);
        if isempty(AllRes)
            fprintf('no units in %s in this file \n',ElLoc);
            return;
        end
        AllRes = round(AllRes*Par.lfpSampleRate/Par.SampleRate)+1;
%         if ~isstr(ElLoc)
%             myi = find(elclu(:,1)==ElLoc(1) & elclu(:,2)==ElLoc(2));
%             AllRes  = AllRes(AllClu==myi);
%             AllClu = ones(length(AllRes),1);
%         end
        sc=1;
        %loop through States
        for s=1:nSt
            %Period = SelectStates(FileBase,States{s},MinPeriod);
            if ~FileExists([FileBase '.sts.' States{s}])
                fprintf('No state %s\n',States{s});
                continue;
            end

            Period = load([FileBase '.sts.' States{s}]);
            PerLen = round(sum(diff(Period,1,2)/Par.lfpSampleRate));
            fprintf('Processing state %s of duration %d seconds\n',States{s},PerLen);

            [Res ind]=SelPerDiscr(AllRes,Period,1,1);
            Clu=AllClu(ind);

            load([FileBase '.thpar.mat'],'ThPh');
            ThPh = SelectPeriods(ThPh,Period,'c',1);
            [dummy PhInd] = histcI(ThPh(Res),linspace(-pi,pi,17));
            phbin = linspace(-pi,pi,17); 
            out(sc).PhBin = (phbin(1:end-1)+phbin(2:end))/2;
            
            eeg = LoadBinary([FileBase '.eeg'],RepCh,Par.nChannels,Period);

            eeg = ButFilter(eeg,2,20/625,'high');
            
            [out(sc).Av out(sc).Std] = TriggeredAv(eeg,Window,Window,Res,Clu);
            for p=1:16
               [out(sc).PhAv(:,:,:,p) out(sc).PhStd(:,:,:,p)] = TriggeredAv(eeg,Window,Window,Res(PhInd==p),Clu(PhInd==p));
            end
            out(sc).nspk = accumarray(Clu,1);
            out(sc).tbin = [-Window:Window]/Par.lfpSampleRate*1000;
            out(sc).fr = FiringRate(Res,Clu,[],1250);
            out(sc).Map = Map;
            out(sc).State = States{s};
            sc= sc+1;
        end %loop through states
        save([FileBase '.' mfilename '.mat'],'out');

        
        
    case 'display'
       load([FileBase '.' mfilename '.mat']);
       tbin = [-25:25]/1.25;
       for i=1:84
           figure(1);clf
           
           mat = NaN*ones(51,96);
           mat(:,RepCh) = out.Av(:,:,i);
           mat = reshape(mat,[51 16 6]);
           mat = reshape(permute(mat,[1 3 2]),[51*6 16]); 
           %PlotTraces(-mat,[],1250,2);axis ij
           imagescnan(mat');
%           PlotTraces96(out.Av(:,:,i),tbin,RepCh,Par,'plot');
           waitforbuttonpress
       end
       
end





