
FileBase = 'sm9603_006-12';
[fMode, ElLoc, FreqRange,States,FrThreshold,Window] = deal('compute','c',[30 150], 'REM',2,2^8);

Par = LoadPar([FileBase '.xml']);
if ischar(States); States = {States}; end;
nSt = length(States);

MinPeriod = 5; %seconds for theta periods selection
if strcmp(ElLoc,'all')
    El = [1:Par.nElecGps];
elseif isstr(ElLoc(1))
    El = find(strcmp(Par.ElecLoc,ElLoc));
else
    El = ElLoc(1);
end

[RepCh ifTetr Info] = RepresentChan(FileBase);
nRepCh = length(RepCh);
eFs = 1250;

MyCh = [3 6 20 33 51 86];
ThCh=  
switch fMode
    case 'compute'
        [AllRes,AllClu,nclu,sph,elclu] = ReadEl4CCG(FileBase,El);
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
        %loop through States
        for s=1:nSt
            Period = SelectStates(FileBase,States{s},MinPeriod);
            if isempty(Period)
                fprintf('No state %s\n',States{s});
                continue;
            end
            PerLen = round(sum(diff(Period,1,2)/1250));
            fprintf('Processing state %s of duration %d seconds\n',States{s},PerLen);

            [Res ind]=SelPerDiscr(AllRes,Period,1,1);
            Clu=AllClu(ind);

            eeg = LoadBinary([FileBase '.eeg'],MyCh,Par.nChannels,Period)';
            
            keyboard
            
            fth =  MTFilter(eeg, [7 2], eFs, 2^10);
            phth = angle(hilbert(fth));
          
            fgam =  MTFilter(eeg, [70 5], eFs, 2^9);
            phgam = angle(hilbert(fgam));
            
        end
        
end
