% function [MySeg, AllStates] = SelectStates(FileBase, StateLabel, MinLen, Display)
% finds the segments that correspond to the staetes in StateLabel (cell array if>1)
% MinLen - in samples , Display - show the summary of each state
function [MySeg, AllStates] = SelectStates(FileBase, varargin)

[StateLabel, MinLen, Display] = DefaultArgs(varargin,{[], 2*1250, 0});
Par = LoadPar([FileBase '.xml']);
nT = FileLength([FileBase '.eeg'])/2/Par.nChannels;

if FileExists([FileBase '.sts.' StateLabel])
    MySeg = load([FileBase '.sts.' StateLabel]);
    AllStates = StateLabel;
    SegDurations = diff(MySeg,1,2);
    MySeg = MySeg(find(SegDurations>MinLen),:);

elseif strcmp(StateLabel,'SWS') 
    if (FileExists([FileBase '.sts.REM']) | FileExists([FileBase '.sts.RUN']) | FileExists([FileBase '.sts.the']))
        %for SWS - if detection is not done - exclude rem and run from the
        %whole file - crude but for now is ok.
        the = [];
        if FileExists([FileBase '.sts.REM'])
            the = [the; load([FileBase '.sts.REM'])];
        end
        if FileExists([FileBase '.sts.RUN'])
            the = [the; load([FileBase '.sts.RUN'])];
        end
        if FileExists([FileBase '.sts.the'])
            the = [the; load([FileBase '.sts.the'])];
        end

        the = sortrows(the,1);
        newper = reshape(the',[],1);
        if newper(1)==1
            newper = newper(2:end);
        else
            newper =[ 1;newper];
        end
        if newper(end)==nT
            newper = newper(1:end-1);
        else
            newper =[newper;nT];
        end
        MySeg = reshape(newper,2,[])';
    else
        MySeg = [1 nT];
    end

    SegDurations = diff(MySeg,1,2);
    MySeg = MySeg(find(SegDurations>MinLen),:);

elseif strcmp(StateLabel,'THE') & (FileExists([FileBase '.sts.REM']) | FileExists([FileBase '.sts.RUN']))
    MySeg = [];
    if FileExists([FileBase '.sts.REM'])
        MySeg = [MySeg; load([FileBase '.sts.REM'])];
    end
    if FileExists([FileBase '.sts.RUN'])
        MySeg = [MySeg; load([FileBase '.sts.RUN'])];
    end
    MySeg = sortrows(MySeg,1);
    SegDurations = diff(MySeg,1,2);
    MySeg = MySeg(find(SegDurations>MinLen),:);
% 
%     
% elseif FileExists([FileBase '.states.res'])
%     if isstr(StateLabel)
%         StateLabel = {StateLabel };
%     end
% 
%     Seg = load([FileBase '.states.res']);
%     SegLen =  load([FileBase '.states.len']);
%     SegClu =  load([FileBase '.states.clu']);
% 
%     Seg2 = Seg(find(SegLen==2),:); Seg2=Seg2(1,:);
%     MinSeg = Seg2(2)-Seg2(1);
%     MinLen = MinLen ./ MinSeg; %now in windows
%     %Labels = {'SWS', 'IMM', 'REM', 'RUN' , 'AWK', 'HVS' , 'ART' ,'DES'};
%     load([FileBase '.statelabels.mat']);
%     MyLabels =[];
%     for i=1:length(StateLabel)
%         thislabel    = find(strcmp(Labels, StateLabel{i}));
%         if ~isempty(thislabel)
%             MyLabels(end+1) = thislabel;
%         end
%     end
%     if isempty(MyLabels)
%         error('state label is wrong');
%         return
%     end
%     MySegInd = find(SegLen>MinLen & ismember(SegClu,MyLabels));
%     if isempty(MySegInd)
%         %fprintf('no some of the states in that file or too short\n');
%         MySeg =[];
%         return
%     end
% 
%     MySeg = Seg(MySegInd,:);
% 
%     Shift = round(MinSeg/2);
%     MySeg = [MySeg(:,1) - Shift MySeg(:,2) + Shift ];
%     AllStates = Labels(MyLabels);
else
    MySeg = [];
   
end