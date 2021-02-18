%
% a function to compute and display the modulation of Ripples by oscill phase
%
% function: OutArgs = PhaseModulation_Gamma(FileBase,  fMode, varargin :: FreqRange, Channel)
%
% fMode: 'compute'; 'display'; 
% FreqRange: [x1 x2] Hz 
% Channel: takes eeg from channel 1; (starts at 0 in the xml!!!)

% The output is a structure OutArgs containing info on phase modulation


function OutArgs = PhaseModulation_Gamma(FileBase, fMode, varargin)

[ FreqRange, Channel] = DefaultArgs(varargin,{[3 6], 15});


switch fMode
    case 'compute'
        
        
        %%% Loading Gamma
        
        load([FileBase '.spw.mat']);
        rips = Rips.t;

        
        %%% Loading the states

        SWS = load([FileBase '.sts.SWS' ]);
       
        
        %%% Loading eeg and getting good segments and the spikes in the segments
        
        Par=LoadPar([FileBase '.xml']);
        eeg  = LoadBinary([FileBase '.eeg'], Channel, Par.nChannels);
        [eeg, orind]= SelectPeriods(eeg(:),SWS,'c',1);
        
       
        
        %%% Filtering the eeg in the Frequency range and getting continuous phase
        
        feeg = ButFilter(eeg, 2, FreqRange/(Par.lfpSampleRate/2),'bandpass');
        hilb = hilbert(feeg);
        ph = angle(hilb);
        
        
        %%% Rayleigh statistic for phase modulation
        
        phstats = PhaseOfTrain(ph(rips),ones(length(rips),1),[],[],1,0,1);
              
        
        %%% Saving the output as OutArgs structure array
        
        OutArgs.cphstats  = CatStruct(phstats);
        OutArgs.ChannelUsed = Channel;
        
            save([FileBase '.' mfilename  '.mat'],'OutArgs');
        

    case 'display'
        

            load([FileBase '.' mfilename  '.mat']);


        figure('name','Phase Modulation - Ripples','NumberTitle','off')
        
        subplot(121)
        [xl, yl]=pol2cart(1, 1); c = compass(xl,yl); set(c, 'Visible', 'off'); hold on % R lim 1
        polar(OutArgs.cphstats.th0, OutArgs.cphstats.R,'o');
        title(['Ripples Phase Modulation']);
        
        subplot(122)
        bar((OutArgs.cphstats.phbin*180/pi)', OutArgs.cphstats.phhist')
        xlim([-200 600]); xlabel(['Phase, deg']);
        title('Ripples Phase Modulation Histogram'); axis tight; hold on
        
end
