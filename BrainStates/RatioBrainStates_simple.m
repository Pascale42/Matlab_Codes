% function RatioBrainStates_simple(FileBase, ch_hpc, ch_cx,)
%
% ch_hpc = hippocampus channel
% ch_cx = cortex channel
%

function RatioBrainStates_simple(FileBase, ch_hpc, ch_cx)


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
powdelta_c=mean(y(:,3:20, 2, 2),2);

 
figure('name', FileBase, 'NumberTitle','off')
subplot 411; imagesc(t,f,20*log10(abs(y(:,:,1,1))+eps)'); axis xy; title('Hpc')
subplot 412; imagesc(t,f,20*log10(abs(y(:,:,2,2))+eps)'); axis xy; title('Cortex')
subplot 413; plot(t, powtheta_h./powdelta_c); axis tight; title ('theta hpc / delta cortex >> REM')
subplot 414; plot(t, powdelta_c./powtheta_h); axis tight; title ('delta cortex / theta hpc >> SWS')
hold on

% kmeans
R=Smooth(powdelta_c./powtheta_h, 6);
plot(t, R, 'LineWidth',2)
idx = kmeans(R, 2);
plot(t, idx*(-10))

[~, on, off]=SchmittTrigger(idx, 1.5, 1.5);

save([FileBase '.' mfilename '.mat'], 'y','f', 't', 'powtheta_h', 'powdelta_c', 'on', 'off')
savefig([FileBase '.' mfilename '.fig']);
