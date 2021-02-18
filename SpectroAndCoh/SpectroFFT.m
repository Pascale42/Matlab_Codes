%
% function Spt = SpectroFFT(FileBase, State, fMode, Channel, ResCoef,
% FreqRange, Overwrite);
%
% fModes: compfreq, dispfreq, comptime, disptime
%

function Spt = SpectroFFT(FileBase, State, fMode, Channel, ResCoef, FreqRange, Overwrite);

switch fMode

  case 'compfreq'

    if exist([FileBase '.' mfilename State '.mat'])>0 & ~Overwrite
      sprintf('Already computed!');return
    end

    % loading;
    Par = LoadPar([FileBase '.xml']);
    STA = load([FileBase '.sts.' State]);


    % loads eeg and resample
    eeg = LoadBinary([FileBase '.eeg'], Channel, Par.nChannels);
    eeg = SelectPeriods(eeg(:),STA,'c',1);
    EEG = resample(eeg,1,ResCoef);

    Spt.FreqRange = FreqRange;
    Spt.reSampleRate = Par.lfpSampleRate/ResCoef;
    Spt.Channel = Channel;


    % for coherence & spectrogram
    if strcmp(State,'THE')==1; WinLengthSec=5; end
    if strcmp(State,'GAMTHE')==1 | strcmp(State,'GAMSWS')==1; WinLengthSec=0.5; end
    if strcmp(State,'SWS')==1; WinLengthSec=15; end

    Spt.WinLengthSample = 2^round(log2(WinLengthSec*Spt.reSampleRate));
    Spt.nFFT = 2*Spt.WinLengthSample;

    % save parameters
    Spt.EEG = EEG;
    Spt.State = State;
    Spt.Par = Par;

    % computes
    Spt.weeg = WhitenSignal(Spt.EEG,Spt.reSampleRate*2000,1);
    [y,f,phi]= mtchd(Spt.weeg,Spt.nFFT,Spt.reSampleRate,Spt.WinLengthSample,[],3,'linear',[],Spt.FreqRange);

    Spt.mtchd.y = y;
    Spt.mtchd.f = f;
    Spt.mtchd.phi = phi;

    save([FileBase '.' mfilename State 'Freq.mat'], 'Spt');
    

    
    
    case 'dispfreq'
    load([FileBase '.' mfilename State 'Freq.mat']);

   
    figure(981160)
          plot(Spt.mtchd.f,20*log10(abs(Spt.mtchd.y)+eps),'r'); grid on;
          ylabel('psd (dB)'); xlabel('Frequency'); grid on; % power spectrum density
    
          
          
    
    
    
    
    
    
    case 'comptime'

    if exist([FileBase '.' mfilename State '.mat'])>0 & ~Overwrite
      sprintf('Already computed!');return
    end

    % loading;
    Par = LoadPar([FileBase '.xml']);
    STA = load([FileBase '.sts.' State]);


    % loads eeg and resample
    eeg = LoadBinary([FileBase '.eeg'], Channel, Par.nChannels);
    eeg = SelectPeriods(eeg(:),STA,'c',1);
    EEG = resample(eeg,1,ResCoef);

    Spt.FreqRange = FreqRange;
    Spt.reSampleRate = Par.lfpSampleRate/ResCoef;
    Spt.Channel = Channel;


    % for coherence & spectrogram
    if strcmp(State,'THE')==1; WinLengthSec=5; end
    if strcmp(State,'GAMTHE')==1 | strcmp(State,'GAMSWS')==1; WinLengthSec=0.5; end
    if strcmp(State,'SWS')==1; WinLengthSec=15; end

    Spt.WinLengthSample = 2^round(log2(WinLengthSec*Spt.reSampleRate));
    Spt.nFFT = 2*Spt.WinLengthSample;

    % save parameters
    Spt.EEG = EEG;
    Spt.State = State;
    Spt.Par = Par;
    Spt.weeg = WhitenSignal(Spt.EEG,Spt.reSampleRate*2000,1);

    % computes
    [y,f,t, phi]= mtchglong(Spt.weeg,Spt.nFFT,Spt.reSampleRate,Spt.WinLengthSample,[],3,'linear',[],Spt.FreqRange);

    Spt.mtchd.y = y;
    Spt.mtchd.f = f;
    Spt.mtchd.t = t;
    Spt.mtchd.phi = phi;

    save([FileBase '.' mfilename State 'Time.mat'], 'Spt');

    
  
          
          
          
          
    
  case 'disptime'
    load([FileBase '.' mfilename State 'Time.mat']);

   
    figure(891160)
    imagesc(Spt.mtchd.t,Spt.mtchd.f,(eps+20*log10(abs(Spt.mtchd.y)))');axis xy; colorbar 
          
          xlabel('Time'); ylabel('Frequency'); 
      
        
end