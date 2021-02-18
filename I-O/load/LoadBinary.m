%function [Data OrigIndex]= LoadBinary(FileName, Channels, nChannels, :: Periods, State, Resample, intype, outtype)
%
% INPUTS::
%                   Channels - list of channels to load starting from 1
%                   nChannels - number of channels in a file, will be read from par/xml file if present
%                   Periods : [beg1 end1; beg2 end2....] (in sampling rate of the file you are loading)
%                   STS : allows to load State (must have a .sts. file)
%                   Resample: resampling factor to get downsampled Data
%                   intype/outtype - data types in the file and in the matrix to load to by default assume input file is eeg/dat = int16 type (short int)
%                                               output is single (to save space) unless it is version 6.5 that cannot handle single type
% OUTPUTS::
%                    Data; each column is a channel
%                   OrigIndex: then returns the original samples index that samples in Data correspond
% to , so that you can use it for future spikes and other point process analysis
%
% By Anton Sirota and Pascale Quilichini


function [data, OrigIndex]= LoadBinary(FileName, Channels, varargin)
if ~FileExists(FileName)
    error('File %s does not exist or cannot be open\n',FileName);
end

lastdot =strfind(FileName,'.');
FileBase=FileName(1:lastdot(end)-1);
if FileExists([FileBase '.xml']) || FileExists([FileBase '.par'])
    Par = LoadPar([FileBase]);
    nChannels = Par.nChannels;
else
    nChannels = 0;
end
[nChannels, Periods, State, Resample, intype, outtype] = DefaultArgs(varargin,{nChannels,[], [],1,'int16','double'});

if ~nChannels
    error('nChannels is not specified and par/xml file is not present')
end

ver = version; ver = str2num(ver(1)); %#ok<ST2NM>
if ver<7;     outtype ='double';end

PrecString =[intype '=>' outtype];
fileinfo = dir(FileName);
% get file size, and calculate the number of samples per channel
nSamples = ceil(fileinfo(1).bytes /datatypesize(intype) / nChannels);

filept = fopen(FileName,'r');

if ~isempty(Periods) && isnumeric(Periods)
    
    if Resample>1
        data = feval( outtype, zeros( length(Channels), sum( ceil((diff(Periods,1,2)+1)/nChannels/Resample) ) ) );
    else
        data = feval( outtype, zeros( length(Channels), sum(diff(Periods,1,2)+1) ) );
    end
end

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

% Old method 2
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
    
    %Pascale 15 Jan 2016
    if ~isempty(State)
        if exist([FileBase '.sts.' State], 'file')
            STA=load([FileBase '.sts.' State]);
            [data, OrigIndex]=SelectPeriods(data(:,:), STA, 'c', 1);
        else
            disp('No sts file for this state, eeg not cut according to STS');
        end
    end
    
end

fclose(filept);



