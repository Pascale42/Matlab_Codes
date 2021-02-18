function Par = LoadEegPar(FileName)
% LoadEegPar(FileName)
% used for eeg analysis 
% loads the specified eeg.par file and returns a structure with these elements:
%   
% .FileName      -> name of file loaded from
% .nElec   -> number of electrodes 
% for Electrode i:
% .ElecChNum(i) -> number of channels that are good in this electrode
% .EleChannels{i}  -> a cell array element giving the channels ID in this electrode
%                    e.g. if .Electrode{3} = [2 3 4 5], electrode 3
%                    is a tetrode for channels 2 3 4 and 5.
% .ElecLoc{i} ->  h/c c for hippocampus/cortex
% .ElecLocType{i} -> 1,3 - CA1/Ca3 in h; number of layer in cx

FileName = [FileName '.eeg.par'];
fp = fopen(FileName, 'r');
Par.FileName = FileName;

Line = fgets(fp);
A = sscanf(Line, '%d');
Par.nElec = A(1);


% read in ElectrodeGroup
for i=1:Par.nElec
    Line = fgets(fp);
    A = sscanf(Line, '%d',1);
    Par.ElecChNum(i) = A(1);
    A = sscanf(Line, '%d');
    Par.ElecChannels{i} = A(2:Par.ElecChNum(i)+1);
    B=sscanf(Line, '%s');
    if (isletter(B(end-1)))
        Par.ElecLoc{i} = B(end-1);
        Par.ElecLocType{i} = B(end);
    else
        Par.ElecLoc{i} = B(end);
        Par.ElecLocType{i} = '';
    end
end;
fclose(fp);
%Par