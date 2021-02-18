function MakeBipolar2(FileBase)
 
Par=LoadPar([FileBase '.xml']);
  
% Loads all channels
EEG=LoadBinary([FileBase '.eeg'], [Par.AnatGrps(1).Channels+1 Par.AnatGrps(2).Channels+1], Par.nChannels);
%EEG=LoadBinary([FileBase '.eeg'], Par.AnatGrps(1).Channels+1, Par.nChannels);
 
% EEGb=zeros(size(EEG,1),size(EEG,2)-4);
EEGb=zeros(size(EEG,1),size(EEG,2)-1); % Arnold EMG Bipol only <<< Pascale

 
% Regular LFP channels
EEGb(:,[1:17])= EEG(:,[1:17]); % Arnold  <<< Pascale

 
% Make Bipolars
% EEGb(:,2) = EEG(:,3) - EEG(:,2); % septum median
% EEGb(:,3) = EEG(:,5) - EEG(:,4); % striatum DL
% EEGb(:,4) = EEG(:,7) - EEG(:,6); % striatum DM
% EEGb(:,7) = EEG(:,11) - EEG(:,10); % Reuniens
 
% EEGb(:,14) = EEG(:,19) - EEG(:,18); % EMG Arnold (Ma?) NIET!
 EEGb(:,18) = EEG(:,19) - EEG(:,18); % Arnold  <<< Pascale; ca peut aussi etre EEGb(:,18) = EEG(:,18) - EEG(:,19); au choix
 
 % Save the new file
 SaveBinary([FileBase '.bip.eeg'], EEGb, 'precision', 'int16', 'mode', 'new')
 
disp('A new EEG BIPOLAR file was created; think about getting the .xml file updated to the new total number of channels!') % <<< C'EST A DIRE UN XML AVEC 18 CANAUX
