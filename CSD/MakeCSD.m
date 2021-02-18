%function csd = MakeCSD(filename,Channels, step, iffilter,Periods, IfIn)
%computes the current source density map on the data in filename
% for the channels 'ch' with 'step' in derivative calculation
% filtering out hightfreq. if iffilter is 1
% Periods and IfIn allows restrict it ot some time intervals only
% the resultand  csd will be a concattenation of the csd in these periods
% Channels = array of size Channels x Locations,
% Channels should contain 1+step*2 consequitive chennel numbers centered
% on channel of the location you want to look at
% so for 3 locations and step 2 you provide 5 by 3
function csdout = MakeCSD(filename,varargin)
[Channels, step,iffilter,Periods, IfIn] =DefaultArgs(varargin,{[],2,1,[], 0});
if isstr(filename)
    par = LoadPar([filename '.par']);
    if isempty(Channels)
        Channels = [1:par.nChannels];
    elseif size(Channels,1)==1 & size(Channels,2)>1
        Channels = Channels';
    end
    nLocations = size(Channels,2);
    nChannels = size(Channels,1);
    LoadChannels = reshape(Channels,nLocations*nChannels,1);
    [LoadChannels dummy Ind] = unique(LoadChannels);
    LoadedChannels = reshape(Ind,nChannels, nLocations);
    eeg = readmulti([filename '.eeg'],par.nChannels,LoadChannels);
    if ~isempty(Periods)
        eeg = SelectedPeriods(eeg,Periods,'c',IfIn);
    end
    if iffilter
        csd = FirFilter(eeg);    
    else 
        csd = eeg;
    end
else
    if size(filename,1)<size(filename,2)
        filename =filename';
    end
    if isempty(Channels)
        Channels = [1:size(filename,2)];
    end
    nLocations = size(Channels,2);
    nChannels = size(Channels,1);
    LoadedChannels = Channels;
    if iffilter
        csd = FirFilter(filename);
    else
        csd = filename;
    end
end
for i=1:nLocations
    nChannels =length(LoadedChannels(:,i));
    csdch=[step+1:nChannels-step];
    csdout(:,i) = -mean(csd(:,csdch+step) - 2*csd(:,csdch) + csd(:,csdch-step),2);
end
if (nargout<1)
    bsave([filename '.csd'], csdout','short');
end
