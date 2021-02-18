 %
% Unit Field Multitaper Coherence Density
% function OutArgs = UF_CoherenceDensity(FileBase, State, fMode, Channel, CM, SaveFig)
%
% State: THE SWS both
% SaveFig: 1 to save report
% fMode: 'compute' and 'display' modes
% CM = gives the subpop of neurons [


function OutArgs = UF_CoherenceDensity(FileBase, State, fMode, Channel, CM, SaveFig);

switch fMode
    
    case 'compute'
        
        % parameters
        Par = LoadPar([FileBase '.xml']);
        if ~exist('CM','var')
            CM=[];
        end
        
        if strcmp('both', State)
            both={'THE';'SWS'};
            for n=1:2
                State=both{n};
                UF_CoherenceDensity_aux(FileBase, State, Channel,CM,Par)
            end
        else
            UF_CoherenceDensity_aux(FileBase, State, Channel,CM,Par)
        end
        
    case 'display'
        if ~exist('CM','var')
            load([FileBase '.' mfilename State '.mat']);
        else
            load([FileBase '.' mfilename State '-cm.mat']);
        end
        
        figure
        imagesc(OutArgs.ufchd.f,[1:length(OutArgs.Map)],squeeze(OutArgs.ufchd.y(:,1,2:end))');
        set(gca,'YTick',[1:size(OutArgs.Map,1)]); set(gca,'YTickLabel',num2str(OutArgs.Map(:,2:3)));
        xlabel('Hz'); title(['Units-Field ' num2str(OutArgs.Channel) ') coherence density for ' FileBase ' - ' num2str(State)]);
        
        if exist('SaveFig','var')
            Par = LoadPar([FileBase '.xml']);
            reportfig(gcf,['UFCoherenceDensity_' State], 0, ['File ' Par.FileName ', State ' State ', ' num2str(OutArgs.FreqRange) ' Hz'],150);
        end
       
end

function UF_CoherenceDensity_aux(FileBase, State, Channel,CM, Par)

STA = load([FileBase '.sts.' State]);
OutArgs.FreqRange=[0.5 20];
ResCoef=1;
OutArgs.reSampleRate = Par.lfpSampleRate/ResCoef;

% loads given channel and resample

if FileExists([FileBase '.eeg.THE.Hpc1Reu.mat'])
    load([FileBase '.eeg.THE.Hpc1Reu.mat']);
    eeg=eeg(:,Channel);
else
    eeg = LoadBinary([FileBase '.eeg'],Channel, Par.nChannels);
    eeg = SelectPeriods(eeg(:),STA,'c',1);
    eeg = resample(eeg,1,ResCoef);
end

% loading Clusters;

[T,G,Map]=LoadCluRes(FileBase);
T = round(T/(Par.SampleRate/Par.lfpSampleRate));
[Tsta Indsta] = SelPerDiscr(T, STA, 1, 1);
Gsta = G(Indsta);
uClu= unique(Gsta);
Map = [uClu Map(uClu,2:3)];
Tsta = round(Tsta/ResCoef);
Tsta(Tsta==0) =1;
Tsta(Tsta>length(eeg)) = length(eeg);
if ~isempty(CM)
    [~,~, idClus] = Intersection(CM, Map(:,2:3));
    Gsta=Gsta(ismember(Gsta,Map(idClus,1))); %!!! ici valeurs 
    Tsta=Tsta(ismember(Gsta,Map(idClus,1)));
    OutArgs.Map=Map(idClus,:); % ici indices
else
    OutArgs.Map=Map;
end
OutArgs.Channel=Channel;

% for coherence & spectrogram
if strcmp(State,'THE')==1; WinLengthSec=5; end
if strcmp(State,'SWS')==1; WinLengthSec=15; end

WinLengthSample = 2^round(log2(WinLengthSec*OutArgs.reSampleRate));
nFFT = 2*WinLengthSample;

[OutArgs.ufchd.y, OutArgs.ufchd.f, OutArgs.ufchd.phi, OutArgs.ufchd.yerr, OutArgs.ufchd.phierr, OutArgs.ufchd.phloc, OutArgs.ufchd.pow]=...
    mtptchd(eeg, Tsta, Gsta, nFFT,OutArgs.reSampleRate,WinLengthSample,[],3,'linear',[],OutArgs.FreqRange);

if ~isempty(CM)
    save([FileBase '.' mfilename State '-cm.mat'], 'OutArgs', 'CM');
else
    save([FileBase '.' mfilename State '.mat'], 'OutArgs', 'CM');
end
