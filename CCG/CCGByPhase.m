%function CCGByPhase(FileBase,fMode,Overwrite, States,El, MaxLag)
% Computes the pairwise ccg for spikes from electrodes El (could be
% string) as a function of theta phase
function out = SpikesCoh(FileBase,varargin)

par = LoadPar([FileBase '.xml']);
[fMode,Overwrite,States,El,MaxLag] = DefaultArgs(varargin,{'compute',1, {'REM','RUN'},'c',20});

if isstr(El)
    El = find(strcmp(par.ElecLoc,El));
end
EegFs =par.lfpSampleRate;
SpikesFs = 1e6/par.SampleTime;

if isstr(States)     States ={States};  end

load([FileBase '.thpar.mat'],'ThPh');
switch fMode
    case 'compute'
        if ~FileExists(SaveFn) | Overwrite
            [Res,Clu,Map] = LoadCluRes(FileBase, El);
            Res = round(Res*EegFs/SpikesFs);
            
            wcnt=1;
            for where = 1:length(States) %loop through states
                Segments = SelectStates(FileBase, States{where},EegFs*2);
                if isempty(Segments)
                    fprintf('No %s periods. Return empty output\n',States{where});
                    continue;
                end
                fprintf('State : %s - PROCESSING\n',States{where});
                %            eeg = LoadBinary([FileBase '.eeg'],Channels, par.nChannels,Segments)';

                %select spikes that are in that State
                [myRes ind] = SelectPeriods(Res,Segments,'d',1,1);
                myClu = Clu(ind);
                [uClu dummy myClu] = unique(myClu);
                myMap =Map(uClu,:);
                uClu = unique(myClu);
                [ccg tbin pairs] = CCG(myRes,myClu,floor(EegFs/1000*MaxLag*2),0,EegFs,uClu,'count');
%                [ccg tbin pairs] = CCG(myRes,myClu,floor(EegFs/1000),MaxLag,EegFs,[1 2],'count');
                if isempty(pairs) continue; end
                %myPairs = pairs(myClu(pairs(:,1))==1 & myClu(pairs(:,2))==2,:);%???
                tLag = diff(myRes(pairs),1,2)/EegFs*1000;
                myPh = ThPh(myRes(pairs(:,1)));
                phBin = linspace(-pi,pi,9);
                tBin = linspace(-MaxLag,MaxLag,20);
                [dumy PhInd] = histcI(myPh,phBin);
                [dumy LagInd] = histcI(tLag,tBin);
                out(wcnt).cnt = accumarray([LagInd PhInd myClu(pairs(:,1)) myClu(pairs(:,2))],1);
                out(wcnt).Map=myMap;
                out(wcnt).State
                out(wcnt).phbin = (phBin(1:end-1)+phBin(2:end))/2;
                out(wcnt).tbin = (tBin(1:end-1)+tBin(2:end))/2;
                
            end
        end
        
        save([FileBase '.' mfilename '.mat'],'out');
end
% for i=1:84
%     for j=i+1:84
%         mycnt = sq(cnt(:,:,i,j));
%         if sum(mycnt(:)>0)<0.3*prod(size(mycnt)) continue; end
%         imagesc(mycnt)
%         waitforbuttonpress
%     end;
% end
% 
% 
% 
% keyboard

