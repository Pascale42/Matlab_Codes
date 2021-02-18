function VarName3 = importmrk_time(filename, startRow, endRow)
%  importmrk_time :: import marker times and makes them in lfpsample time
%IMPORTFILE Import numeric data from a text file as column vectors.
%   VARNAME3 = IMPORTFILE(FILENAME) Reads data from text file FILENAME for
%   the default selection.
%
%   VARNAME3 = IMPORTFILE(FILENAME, STARTROW, ENDROW) Reads data from rows
%   STARTROW through ENDROW of text file FILENAME.
%
% Example:
%   VarName3 = importfile('20131029_HpcReu.gamma.SWS1.mrk',1, 727);
%

%% Initialize variables.
delimiter = '\t';
if nargin<=2
    startRow = 1;
    endRow = inf;
end

%% Format string for each line of text:
%   column3: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%*s%*s%f%*s%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to format string.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
% % % dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'HeaderLines', startRow(1)-1, 'ReturnOnError', false);
% % % for block=2:length(startRow)
% % %     frewind(fileID);
% % %     dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'HeaderLines', startRow(block)-1, 'ReturnOnError', false);
% % %     dataArray{1} = [dataArray{1};dataArrayBlock{1}];
% % % end
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'HeaderLines', startRow(1)-1, 'ReturnOnError', false);
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'HeaderLines', startRow(block)-1, 'ReturnOnError', false);
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

%% Close the text file.
fclose(fileID);

%% Allocate imported array to column variable names
VarName3 = dataArray{:, 1};
VarName3 =VarName3 (2:end);
VarName3=round(VarName3*1250);

VarName2 = dataArray{:, 2};
VarName2 =VarName2 (2:end);
% VarName3=[VarName2 VarName3];
