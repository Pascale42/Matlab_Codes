% function ToAnywave(FileBase, channels, IR, varargin :: State, labels)
%
% channels to load: indices from 1
% IR: input range in mV; one number if same for all the channels, else provide as many values as channels
% State: if any, string of the state for loading FileBase.sts.STATE file and cutting eeg accordingly
% Labels: {'lfp'} by default in not given otherwise
%
% Output: FileBase.AW.dat et FileBase.AW.ades files for AnyWave

function ToAnywave(FileBase, channels, IR, varargin)

[State, labels] = DefaultArgs(varargin,{[],[]});

Par=LoadPar([FileBase '.xml']);
eeg=GetEegState(FileBase, channels, 0, State);

nch=length(channels);


if length(channels) == 1
    eeg=GetmV(eeg, IR);
    
else
    if length(channels) ~= length(IR)
        IR=repmat(IR, length(channels),1);
    end
    
    for n=1:length(channels)
        eeg(:,n)=GetmV(eeg(:,n), IR(n));
    end
    
end

if isempty(labels)
    labels = repmat({'lfp'}, [length(channels) 1]);
end

mat2ades(eeg',[FileBase '.AW'],1250,labels,'EEG'); 


%% DEPENDENCIES

function varargout = DefaultArgs(Args, DefArgs)
% auxillary function to replace argument check in the beginning and def. args assigment
% sets the absent or empty values of the Args (cell array, usually varargin)
% to their default values from the cell array DefArgs.
% Output should contain the actuall names of arguments that you use in the function

if isempty(Args)
    Args ={[]};
end

if ~iscell(DefArgs)
    DefArgs = {DefArgs};
end
nDefArgs = length(DefArgs);
nInArgs = length(Args);

if (nargout~=nDefArgs)
    error('number of defaults is different from assigned');
end

for n=1:nDefArgs
    
    if (n>nInArgs || isempty(Args{n}))
        varargout(n) = {DefArgs{n}}; %#ok<*CCAT1>
    else
        varargout(n) = {Args{n}};
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%
function [eeg, orind]=GetEegState(FileBase, Channels, save, State)


Par = LoadPar([FileBase '.xml']);

disp('... loading eeg ...')

eeg  = LoadBinary([FileBase '.eeg'], Channels , Par.nChannels);
orind=[];

if ~isempty(State)
    STA = load([FileBase '.sts.' State]);
    [eeg, orind] = SelectPeriods(eeg(:,:),STA,'c',1);
end

if save == 1
    save([FileBase '.eeg.' State '.mat'], 'eeg', 'orind', 'Channels', 'State', 'STA', '-v7.3')
end

%%%%%%%%%%%%%%%%%%%%%%%%%
function [y, ind] = SelectPeriods(x,Periods,varargin)

% selects from the signal (time series  - 'd', continuous - 'c')
% values in/out the Periods ranges (should be the same sampl.rate
% WhereFlag defines in(1) or out(0) (default - 1 = in)
% ifSquash is 1 if you apply it for discrete signal and want to squash the
% gaps to match the sample indexes of continuous signal
%   output: y - new signal, ind - indexes of the old signal

[SigType,WhereFlag, ifSquash] = DefaultArgs(varargin,{'c',1,0});

if isstr(Periods) %#ok<*DISSTR>
    Periods = load(Periods);
end
if isempty(Periods)
    y = x;
    ind = [1:length(x)];
    return
end


Periods(find(Periods(:)==0))=1;

nPeriods  = size(Periods,1);

if ~iscell(x)
    nChannels = min(size(x));
    if size(x,1)==nChannels
        x=x';
    end
    nTimeBins = max(size(x));
    
    if (nargin<3 || isempty(SigType))
        if (nTimeBins < max(Periods(:)))
            SigType = 'd';
        else
            SigType = 'c';
            if (nargout >2)
                error('too many output parameters');
                exit  %#ok<*UNRCH>
            end
        end
    end
    if (SigType =='d')
        Channels = size(x,1);
    end
else
    celli=size(x,1);
    cellj=size(x,2);
end


if (nargin <4 || isempty(WhereFlag) )
    WhereFlag =1;
end



if (SigType == 'd')
    if ~iscell(x)
        [y, ind] = SelPerDiscr(x,Periods,WhereFlag,ifSquash);
    else
        y = cell(celli,cellj);
        ind = cell(celli,cellj);
        for ii=1:celli
            for jj=1:cellj
                [y{ii,jj}, ind{ii,jj}]= SelPerDiscr(x{ii,jj},Periods,WhereFlag,ifSquash);
            end
        end
    end
    
    
end


if (SigType == 'c')
    y=[];
    ind=[];
    
    
    if WhereFlag
        for p=1:nPeriods
            y = [y; x(Periods(p,1):Periods(p,2),:)];
            ind = [ind; [Periods(p,1):Periods(p,2)]'];
        end
    else
        if Periods(1,1)>2
            y = [y; x(1:Periods(1,1)-1,:)];
            ind = [ind; [1:Periods(1,1)-1]'];
        end
        for p=1:nPeriods-1
            y = [y; x(Periods(p,2)+1:Periods(p+1,1)-1,:)];
            ind = [ind; [Periods(p,2)+1:Periods(p+1,1)-1]'];
        end
        
        y = [y; x(Periods(nPeriods,2)+1:end,:)];
        ind = [ind; [Periods(nPeriods,2)+1:size(x,1)]'];
        
    end
    
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%
function [data, OrigIndex]= LoadBinary(FileName, Channels, nChannels, varargin )
%   Channels - list of channels to load starting from 1
%   nChannels - number of channels in a file
%   intype/outtype - data types in the file and in the matrix to load to by default assume input file is eeg/dat = int16 type (short int), 
%   and output is single (to save space) unless it is version 6.5 that cannot  handle single type
%   Works with buffers-works even for huge files
%  OrigIndex then returns the original samples index that samples in Data correspond to, so that you can use it for future spikes and other point process analysis


% Prerequisites

if ~FileExists(FileName)
    error('File %s does not exist or cannot be open\n',FileName);
end

lastdot =strfind(FileName,'.');
FileBase=FileName(1:lastdot(end)-1);

[Periods, intype, outtype,Resample] = DefaultArgs(varargin,{[],'int16','double',1});

if ~nChannels
    error('nChannels is not specified')
end

ver = version; ver = str2num(ver(1)); %#ok<ST2NM>
if ver<7;     outtype ='double';end

PrecString =[intype '=>' outtype];
fileinfo = dir(FileName);


% Get file size, and calculate the number of samples per channel

nSamples = ceil(fileinfo(1).bytes /datatypesize(intype) / nChannels);

filept = fopen(FileName,'r');

if ~isempty(Periods) && isnumeric(Periods)
    

% Resampling
    
    if Resample>1
        data = feval( outtype, zeros( length(Channels), sum( ceil((diff(Periods,1,2)+1)/nChannels/Resample) ) ) );
    else
        data = feval( outtype, zeros( length(Channels), sum(diff(Periods,1,2)+1) ) );
    end
end


% Periods

if isempty(Periods)
    Periods = [1 nSamples];
    if Resample>1
        data = feval(outtype,zeros(length(Channels), ceil(nSamples/Resample)));
    else
        data = feval(outtype,zeros(length(Channels), nSamples));
    end
end

if ~isempty(Periods) && ischar(Periods)
    State=Periods;
    Periods = [1 nSamples];
    if Resample>1
        data = feval(outtype,zeros(length(Channels), ceil(nSamples/Resample)));
    else
        data = feval(outtype,zeros(length(Channels), nSamples));
    end
end


% Original Indexing (OrigIndex)

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


%%%%%%%%%%%%%%%%%%%%%%%%%
function exists = FileExists(FileName)
% FileExists(FileName)
%
% returns 1 if the file exists, 0 otherwise
fp = fopen(FileName, 'r');

if (fp==-1)
	exists = 0;
else
	fclose(fp);
	exists = 1;
end;



%%%%%%%%%%%%%%%%%%%%%%%%%
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
