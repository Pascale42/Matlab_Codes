
function [T, G, Map] = LoadCluResKS(varargin)
% function [T, G, Map] = LoadCluResKS( varargin: myKsDir, SamplingRate, Save)
%
% Get spike times, cluster and Map
%
% outputs:
% - T is length nSpikes vector with time in sample of every spike
% - G is corresponding cluster ID
% to the position of the template it was extracted with)
% - Map is the matrix giving the AnatGroup (row 1),  cluster ID (row 2) and
% on which Channel (the real channel, from 1) its amplitude is the biggest (row 3)
%
% inputs:
% - myKsDir; path to find the data to load; pwd by default
% - SamplingRate; if not parsed, gets it from the xml
% - Save; 0 for not saving, 1 saves as FileBase.CluRes.mat
%
% only for clusters tagged as Good
%

[myKsDir, SamplingRate, Save] = DefaultArgs(varargin,{pwd, [], 0});

% Gets the SamplingRate
if isempty(SamplingRate)
    a=dir;
    for n=1:size(a); if ~isempty(strfind(a(n).name, 'xml')) ; FileName = a(n).name; end; end
    Par = LoadPar(FileName);
    SamplingRate = Par.SampleRate; clear  dir a
end

sp = loadKSdir(myKsDir);
cd(myKsDir)

% get only Good clusters ()
Good = ismember(sp.clu, sp.cids(sp.cgs == 2));
G = sp.clu(Good);
G=double(G); % Gracias Ana!
T = sp.st(Good)*SamplingRate; % go in samples



% FROM templatePositionsAmplitudes.m_____________

% unwhiten all the templates
tempsUnW = zeros(size(sp.temps));
for t = 1:size(sp.temps,1)
    tempsUnW(t,:,:) = squeeze(sp.temps(t,:,:))*sp.winv;
end

% The amplitude on each channel is the positive peak minus the negative
tempChanAmps = squeeze(max(tempsUnW,[],2))-squeeze(min(tempsUnW,[],2));

% The template amplitude is the amplitude of its largest channel (but see
% below for true tempAmps)
tempAmpsUnscaled = max(tempChanAmps,[],2);

% need to zero-out the potentially-many low values on distant channels ...
threshVals = tempAmpsUnscaled*0.3;
tempChanAmps(bsxfun(@lt, tempChanAmps, threshVals)) = 0;
%
%             % ... in order to compute the depth as a center of mass
%             templateDepths = sum(bsxfun(@times,tempChanAmps,sp.ycoords'),2)./sum(tempChanAmps,2);

% __________________________________________________

% channel for each template
[~, max_site]=max(tempChanAmps,[],2);

spikeMaxChan=zeros(length(sp.spikeTemplates), 1);
for n=1:length(sp.spikeTemplates)
    if ~sp.spikeTemplates(n) == 0
        spikeMaxChan(n) = max_site(sp.spikeTemplates(n));
    end
end

% Get only for Good clusters
spikeMaxChan = spikeMaxChan(Good);
spikeTemplates= sp.spikeTemplates(Good);

% Make Map
Map = NaN(length(unique(G)),3);
Map(:,2) = unique(G);

load('chanMap.mat', 'kcoords')
ChanList=[];
for n=1:size(Par.AnatGrps,2)
    ChanList=[ChanList Par.AnatGrps(n).Channels+1];
end


for n=1:length(unique(G))
    current_temp = mode(spikeTemplates(G==Map(n,2)));
    Map(n,3)= ChanList(max_site(current_temp+1)); % channel with max amplitude of template
    chan= max_site(current_temp+1);
    Map(n,1) = kcoords(chan) ; % make AnatGroup
    %     while chan - 32 > 0
    %         Map(n,1) = Map(n,1)+1; chan=chan-32;
    %     end
end

[~,idx] = sort(Map(:,1));
Map = (Map(idx,:));



if Save
    save([FileName(1:end-4) '.CluRes.mat'], 'T', 'G', 'Map')
end