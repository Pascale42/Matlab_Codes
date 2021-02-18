% Function Data = LoadFilEeg(FileBase, State, Channels, tag, Region, FreqRange, ResCoef);
%
% State: THE, SWS
% Region: hpc, ec
%
% saves as: FileBase.LoadFilEeg-Region-State-.mat


function Data = LoadFilEeg(FileBase, State, Channels, tag, Region, FreqRange, ResCoef);


% Load data and parameters file and Resample eeg;

STS = load([FileBase '.sts.' State]);
Par = LoadPar([FileBase '.xml']);

EEG=[]; OrInd=[];

        for n=Channels
            [eeg orind] = LoadBinary([FileBase '.eeg'], n, Par.nChannels,STS);
            eeg = resample(eeg',1,ResCoef); % transpose eeg because ButFilter assumes that eeg vect is a column;
            orind = resample(orind,1,ResCoef);
            EEG=[EEG eeg]; OrInd=[OrInd orind];
        end
        
 % Filtering;

reSampleRate = Par.lfpSampleRate/ResCoef;
feeg = ButFilter(EEG, 2, FreqRange/(reSampleRate/2),'bandpass');


Para.Channels = Channels;
Para.reSampleRate = Par.lfpSampleRate/ResCoef;
Para.ResCoef = ResCoef;
Para.FreqRange = FreqRange;


save([FileBase '.LoadFilEeg-' Region '_' tag '-' State '.mat'], 'EEG', 'Para', 'feeg', 'OrInd')



% SelectPeriods brings eeg into a column!! no need to transpose
