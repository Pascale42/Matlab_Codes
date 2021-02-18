%
% :: USAGE ::
% function data= LoadBinaryData(FileName, Channels, nChannels)
%
% :: INPUTS ::
%   Channels - list of channels to load starting from 1
%   nChannels - number of channels in a file
%
% :: ADDITIONNAL INPUT ARGUMENTS ::
%  intype/outtype - data types in the file and in the matrix to load 
%  by default assume input file is eeg/dat = int16 type (short int), and output is single 
% (to save space) unless it is version 6.5 that cannot handle single type
%  Resample - resampling coefficient for heavy data
% 
% :: ADDITIONNAL OUTPUT ::
% OrigIndex -  returns the original samples index that samples in data correspond to, 
% so that you can use it for future spikes and other point process analysis
%
% Pascale Quilichini; INS - INSERM U1106 - (pascale.quilichini@univ-amu.fr) - 2015
%


function [data, OrigIndex]= LoadBinaryData(FileName, Channels, nChannels, varargin)

% Prerequisite
if ~FileExists(FileName)
    error('File %s does not exist or cannot be open\n',FileName);
end

lastdot =strfind(FileName,'.');
FileBase=FileName(1:lastdot(end)-1);

[intype, outtype,Resample] = DefaultArgs(varargin,{'int16','double',1});

ver = version; ver = str2num(ver(1)); %#ok<ST2NM>
if ver<7;     outtype ='double';end

PrecString =[intype '=>' outtype];
fileinfo = dir(FileName);


% Get file size, and calculate the number of samples per channel
nSamples = ceil(fileinfo(1).bytes /datatypesize(intype) / nChannels);

filept = fopen(FileName,'r');

Periods = [1 nSamples];
data = feval(outtype,zeros(length(Channels), nSamples));


% Buffered way
OrigIndex = [];
nPeriods = size(Periods,1);
buffersize = 400000;
if Resample>1 ; buffersize = round(buffersize/Resample)*Resample; end
totel=0;

for ii=1:nPeriods
    numel=0;
    numelm=0;
    Position = (Periods(ii,1)-1)*nChannels*datatypesize(intype);
    ReadSamples = diff(Periods(ii,:))+1;
    fseek(filept, Position, 'bof');
    
    while numel<ReadSamples
        if numel==ReadSamples; break; end
        [tmp,count] = fread(filept,[nChannels,min(buffersize,ReadSamples-numel)],PrecString);
        data(:,totel+1:totel+ceil(count/nChannels/Resample)) = tmp(Channels,1:Resample:end);
        numel = numel+count/nChannels;
        totel = totel+ceil(count/nChannels/Resample);
    end
    
    OrigIndex = [OrigIndex; Periods(ii,1)+[0:Resample:ReadSamples-1]'];
    data=data';
    
end

fclose(filept);

%% Auxiliary functions

function varargout = DefaultArgs(Args, DefArgs)
% auxillary function to replace argument check in the beginning and def. args assigment
% sets the absent or empty values of the Args (cell array, usually varargin)
% to their default values from the cell array DefArgs.
% Output should contain the actuall names of arguments that you use in the function

% e.g. : in function MyFunction(somearguments , varargin)
% calling [SampleRate, BinSize] = DefaultArgs(varargin, {20000, 20});
% will assign the defualt values to SampleRate and BinSize arguments if they
% are empty or absent in the varargin cell list
% (not passed to a function or passed empty)
if isempty(Args)
    Args ={[]};
end

% if iscell(Args) & isstr(Args{1}) & length(Args)==1
%     Args = Args{1};
% end

if ~iscell(DefArgs)
    DefArgs = {DefArgs};
end
nDefArgs = length(DefArgs);
nInArgs = length(Args);
%out = cell(nDefArgs,1);
if (nargout~=nDefArgs)
    error('number of defaults is different from assigned');
    %keyboard
end
for i=1:nDefArgs
    
    if (i>nInArgs | isempty(Args{i}))
        varargout(i) = {DefArgs{i}};
    else
        varargout(i) = {Args{i}};
    end
end

function exists = FileExists(FileName)
% returns 1 if the file exists, 0 otherwise
fp = fopen(FileName, 'r');

if (fp==-1)
	exists = 0;
else
	fclose(fp);
	exists = 1;
end;

function TypeSize = datatypesize(DataType)
% returns the number of byte for a give format
if strcmp(DataType, 'uint8') | strcmp(DataType, 'int8')
    TypeSize = 8;
elseif strcmp(DataType, 'uint16') | strcmp(DataType, 'int16') | strcmp(DataType, 'short')
    TypeSize = 16;
elseif  strcmp(DataType, 'uint32') | strcmp(DataType, 'int32') | strcmp(DataType, 'long')
    TypeSize = 32;
end
TypeSize = TypeSize/8;
