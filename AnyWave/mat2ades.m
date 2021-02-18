function mat2ades(data,fileName,FS,labels,labelType) 
 
% write in the current folder ADES and DAT files from matrix in MATLAB workspace
% data = matrix of data (nbchannel * time points) - the data have to be in
% microVolt
% fileName = string of the output files without extension ; the ades and
% dat files will have the same name
% FS = sampling rate 
% labels = cell array with channel labels
% labelType : 'EEG' or 'MEG'
% Sophie Chen - January 2014
 
%% generate the ADES file
adesFile = [fileName '.ades'];
 
fid = fopen(adesFile,'wt');
 
fprintf(fid,'%s\r\n','#ADES header file ');
fprintf(fid,'%s','samplingRate = ');
fprintf(fid,'%d\r\n',FS);
fprintf(fid,'%s','numberOfSamples = ');
fprintf(fid,'%d\r\n',size(data,2));
 
for lab = 1:length(labels)
    fprintf(fid,'%s\r\n',[labels{lab} ' = ' labelType]);
end
 
fclose(fid);
 
%% generate the DAT file
 
datFile = [fileName '.dat'];
 
fad = fopen(datFile,'wb');
fwrite(fad, data, 'float32');
fclose(fad);
