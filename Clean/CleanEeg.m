% function CleanEeg(FileBase)
%
% Must create a .sts.Clean file before

% Made in Maeva


function CleanEeg(FileBase)

if ~FileExists([FileBase '.sts.Clean'])
    disp('Run CheckEegStates and Create FileBase ".sts.Clean" first!!!')
    return
end


% load the Par and  the eeg file (all the channels; here isan example for a total of 2 AnatGroups)
 Par=LoadPar([FileBase '.xml']);
 
 disp('Loading and cutting...');
%  Channels=[]; 
% for n=1:size(Par.AnatGrps,2) NON!!!! respecter ordre original des canaux
%   Channels=[Channels Par.AnatGrps(n).Channels+1];
%  end
% 
%  EEG = LoadBinary([FileBase '.eeg'], Channels, Par.nChannels);
 
EEG = LoadBinary([FileBase '.eeg'],  1:Par.nChannels, Par.nChannels);

 % load the Clean segments
 Clean=load([FileBase '.sts.Clean']);
 
 % make dirty segments
 if size(Clean,1) >1
     dirty=zeros(size(Clean,1)+1, 2);
     dirty(1,:) = [1 Clean(1,1)];

     for n=1:size(Clean,1)-1
         dirty(n+1,:)= [Clean(n,2) Clean(n+1,1)];
     end
     dirty(size(Clean,1)+1,:)= [Clean(size(Clean,1),2) size(EEG,1)];
     
 else
     dirty=zeros(2, 2);
     dirty(1,:) = [1 Clean(1,1)];
     dirty(2,:) = [Clean(1,2) size(EEG,1)];   
 end

 % replace eeg dirty parts
 disp('replace eeg dirty parts')
 for n=1:size(dirty,1)
    EEG(dirty(n,1) : dirty(n,2),:) = rand(size(EEG(dirty(n,1) : dirty(n,2),:),1),  size(EEG,2)) .*10^4;
 end

 % Rename old EEG for saving
 disp('Rename old EEG for saving')
copyfile([FileBase '.eeg'], [FileBase '.ini.eeg']);
 
 % Save the new file
 disp('Save the new file')
 SaveBinary([FileBase '.eeg'], EEG, 'precision', 'int16', 'mode', 'new')

 disp('A new EEG cleaned file was saved')

 
 %%  OLD function CleanEeg(FileBase)
% % %
% % % Must create a .sts.Clean file before
% % % Assumes that there are 2 Anat Groups!!
% % %Made in Maeva
% % 
% % 
% % function CleanEeg(FileBase)
% % 
% % if ~FileExists([FileBase '.sts.Clean'])
% %     disp('Run CheckEegStates and Create FileBase ".sts.Clean" first!!!')
% %     return
% % end
% % 
% % 
% % % load the Par and  the eeg file (all the channels; here isan example for a total of 2 AnatGroups)
% %  Par=LoadPar([FileBase '.xml']);
% %  
% %  disp('Loading and cutting...');
% %  EEG=LoadBinary([FileBase '.eeg'], [Par.AnatGrps(1).Channels+1 Par.AnatGrps(2).Channels+1], Par.nChannels);
% %  
% %  % load the Clean segments
% %  Clean=load([FileBase '.sts.Clean']);
% %  
% %  % cut out the Clean segments out of the eegs
% %  EEG=SelectPeriods(EEG(:,:), Clean, 'c', 1);
% %  
% %  % Rename old EEG for saving
% % copyfile([FileBase '.eeg'], [FileBase '.ini.eeg']);
% %  
% %  % Save the new file
% %  SaveBinary([FileBase '.eeg'], EEG, 'precision', 'int16', 'mode', 'new')
% % % %  
% % % % %  % Make a file to cut videos with cleanvid command in bash
% % % % clean=Clean/Par.lfpSampleRate; %  
% % % % for n=1:size(clean,1) 
% % % % clean(n,2)=clean(n,2) - clean(n,1); 
% % % % end %
% % % % 
% % % % fid=fopen([FileBase '.clean.txt'], 'w'); 
% % % % for n = 1:size(clean,1) %
% % % % fprintf(fid,'%f\t',clean(n,:)); fprintf(fid,'\n');
% % % %  end %
% % % % fclose(fid);
% %  
% %  disp('A new EEG Clean file was saved')
