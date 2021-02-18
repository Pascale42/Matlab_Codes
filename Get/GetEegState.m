% [eeg, orind] = GetEegState(FileBase, Channels, Save, State);
%
% State is an option, use [] instead for the whole eeg

function [eeg, orind]=GetEegState(FileBase, Channels, Save, State)


Par = LoadPar([FileBase '.xml']);

disp('... loading eeg ...')

eeg  = LoadBinary([FileBase '.eeg'], Channels , Par.nChannels);
orind=[];

if ~isempty(State)
    STA = load([FileBase '.sts.' State]);
    [eeg, orind] = SelectPeriods(eeg(:,:),STA,'c',1);
end

if Save == 1
    save([FileBase '.eeg.' State '.mat'],'-v7.3', 'eeg', 'orind', 'Channels', 'State', 'STA')
end



