function CleanDat(FileBase)

% %En amont:
%
% % Do the CheckEegStates on your original eeg file
% CheckEegStates(FileBase,'ART',FreqRange,ChannelsChosen);

% A faire tourner sur le Cluster >>>


% load the Par
Par=LoadPar([FileBase '.xml']);


% load the Artefact segments
dirty=load([FileBase '.sts.ART']);

dat=LoadBinary([FileBase '.dat'], 1:Par.nChannels, Par.nChannels);

% for each channel
for ch=1:Par.nChannels
    
    for n=1:size(dirty,1)
        dat(dirty(n,1) : dirty(n,2),ch) = rand(length(dat(dirty(n,1) : dirty(n,2),:)),1) .* std(dat(:,ch)); 
    end
    
end


% Save the new file
SaveBinary([FileBase '_cleaned.dat'], dat,  'precision', 'int16', 'mode', 'new')