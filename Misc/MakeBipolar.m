function MakeBipolar(FileBase)

Par=LoadPar([FileBase '.xml']);
  
% Loads all channels
EEG=LoadBinary([FileBase '.eeg'], [Par.AnatGrps(1).Channels+1 Par.AnatGrps(2).Channels+1],  Par.nChannels);

EEGb=zeros(size(EEG,1),size(EEG,2)-4);

% Regular LFP channels
%EEGb(:,[1 5 6 8:17])=EEG(:,[1 8 9 12:21]);

EEGb(:,[1 5 6 8:13])=EEG(:,[1 8 9 12:17]); % Arnold (Maë)

% Make Bipolars
EEGb(:,2) = EEG(:,3) - EEG(:,2); % septum median
EEGb(:,3) = EEG(:,5) - EEG(:,4); % striatum DL
EEGb(:,4) = EEG(:,7) - EEG(:,6); % striatum DM
EEGb(:,7) = EEG(:,11) - EEG(:,10); % Reuniens

EEGb(:,14) = EEG(:,19) - EEG(:,18); % EMG Arnold (Maë)


 % Save the new file
 SaveBinary([FileBase '.bip.eeg'], EEGb, 'precision', 'int16', 'mode', 'new')
 
 disp('A new EEG BIPOLAR file was created; think about getting the .xml file!')