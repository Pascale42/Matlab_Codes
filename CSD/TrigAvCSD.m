%function TrigAvCSD(filename, Triggers,  win, WhichCh, Periods, IfIn, Fs, step, nChannels, ColorRange,method,orient)
% Triggers in sampling rate of the file to trigger average,
% win - window in msec
% method : 1, to load all channels for each segment around trigger (if nChannels/nTriggers is large)
%               2, to load each channel and trigger average it separately (use that if nChannels/nTriggers is small)
% orient - 'v' if vertical, and 'h' for horizontal if only one of the
% WhichCh or Triggers is a cell array, otherwise igonred
function [TrAv, t, TrStd] =TrigAvCSD(filename, Triggers, varargin)

[win, WhichCh, Periods, IfIn, Fs, step, nChannels, ColorRange,method, orient]  = DefaultArgs(varargin, { 200, [], [], 1, 1250, 2, [], [],2, 'v' });
if iscell(WhichCh)
    nShanks = length(WhichCh);
else
    WhichCh = {WhichCh};
    nShanks =1;
end
if iscell(Triggers)
    %get rid of empty triggers (may come in batch mode)
    newT ={};
    for i=1:length(Triggers)
        if ~isempty(Triggers{i})
            newT{end+1} = Triggers{i};
        end
    end
    Triggers = newT;
    nTriggers = length(Triggers);
else
    nTriggers = 1;
    Triggers = {Triggers};
end
nTriggers = length(Triggers);
if isstr(Periods)
    Periods = load(Periods);
    Periods = round(Periods/16);
end

%exclude points from outside Periods if needed
if ~isempty(Periods)
    newTr ={};
    for tr=1:nTriggers
        if ~isempty(Periods)
            Tr = SelectPeriods(Triggers{tr}, Periods, 'd', IfIn);
        else
            Tr = Triggers{tr};
        end
        if ~isempty(Tr)
            newTr{end+1} = Tr;
        end
    end
    nTriggers = length(newTr);
    Triggers = newTr;
end

for sh=1:nShanks
    for tr=1:nTriggers
        [TrAv{sh, tr}, TrStd{sh, tr}, t]  = TriggeredAvM(filename, Triggers{tr}, win, Fs, nChannels,method);
        TrAv{sh, tr} = TrAv{sh, tr}(:,WhichCh{sh});
    end
end

if nTriggers>1 & nShanks>1 
    orient = [];
end
%reorder the TrAv according to the WhichCh order
if size(WhichCh,1)==1
    NewTrAv ={};
    for sh=1:nShanks
        for tr=1:nTriggers
            NewTrAv{tr,sh} = TrAv{sh, tr};
        end
    end
    TrAv = NewTrAv;
end

if nargout<1
    %  clf
    PlotCSD(TrAv, t, [], step, ColorRange,1,1,[],[],orient);
end
