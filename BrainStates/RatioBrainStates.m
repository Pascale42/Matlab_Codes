% function RatioBrainStates(FileBase, ch_hpc, ch_cx, emg)
%
% ch_hpc = hippocampus channel
% ch_cx = cortex channel
% emg = emg channel OR give the vector of lfp emg
%

function RatioBrainStates(FileBase, ch_hpc, ch_cx, emg)


% Loading
Par=LoadPar([FileBase '.xml']);
eeg=LoadBinary([FileBase '.eeg'], [ch_hpc ch_cx] , Par.nChannels);

% Pre-whitening
 eeg = WhitenSignal(eeg,Par.lfpSampleRate*2000, 2);

 % Parameters for TF
 WinLengthSec = 4;
 WinLengthSample = 2^round(log2(WinLengthSec*Par.lfpSampleRate));
 nFFT = 2*WinLengthSample;

 [y, f, t] = mtchglong(eeg(:,1:2),nFFT,Par.lfpSampleRate,WinLengthSample,[],3,'linear',[],[0.1 70]);

% delta = 0.5-1.7 Hz > f(3:11)
% theta = 5.5-8.5 Hz > f(36:56)
% gamma = 52-65 Hz > f(341:426)
% beta = 20-25 Hz f(132:163)

% canal hpc
powtheta_h=mean(y(:,36:56, 1, 1), 2);

% canal cortex
powdelta_c=mean(y(:,3:11, 2, 2),2);

% canal emg
if length(emg) == 1
eeg_emg=LoadBinary([FileBase '.eeg'], emg , Par.nChannels);
else
    eeg_emg=emg; clear emg
end
 [yemg, femg, temg] = mtchglong(eeg_emg,nFFT,Par.lfpSampleRate,WinLengthSample,[],3,'linear',[],[150 250]);
 pow_emg=mean(yemg, 2);

 
figure;
subplot 611; imagesc(t,f,20*log10(abs(y(:,:,1,1))+eps)'); axis xy; title('Hpc')
subplot 612; imagesc(t,f,20*log10(abs(y(:,:,2,2))+eps)'); axis xy; title('Cortex')
subplot 613; imagesc(temg,femg,20*log10(abs(yemg(:,:))+eps)'); axis xy; title('EMG')
subplot 614; plot(t, powtheta_h./powdelta_c); axis tight; title ('theta hpc / delta cortex >> REM')
subplot 615;plot(t, powtheta_h.*pow_emg); axis tight; title ('theta hpc * EMG >> THE RUN')
subplot 616; plot(t, powdelta_c./pow_emg); axis tight; title ('delta cortex / EMG >> SWS')


save([FileBase '.' mfilename '.mat'], 'y','f', 't', 'pow_emg', 'pow_theta_h', 'powdelta_c', 'temg', 'yemg', 'femg')
savefig([FileBase '.' mfilename '.fig']);