
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% PARAMETERS %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



WinLengthSec=5; % for theta

WinLengthSec=0.5; % for gamma and up

WinLengthSec=15; % for 1 Hz
       
WinLengthSample = 2^round(log2(SR));
nFFT = 2*OutArgs.WinLengthSample;


%%%%%%%%%%%% Whiten Signal %%%%%%%%%%%%

weeg = WhitenSignal(eeg,SR*2000,2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% COHERENCE DENSITY %%%%%%%%%%%%     %%%%%%%%%%%% MTCHD %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


      %%%%%%%%%%%% compute %%%%%%%%%%%%

[y,f,phi]= mtchd(eeg,nFFT,SR,WinLengthSample,[],3,'linear',[],FreqRange);

[y, f, phi, yerr, phierr, phloc, pow]=...
            mtptchd(eeg, [], [], nFFT,SR,WinLengthSample,[],3,'linear',[],FreqRange)

        %%%%%%%%%%%% plot %%%%%%%%%%%%      


         figure(24)
        for Ch1=1:nChannels
            for Ch2 = Ch1:nChannels
                subplot(nChannels, nChannels, Ch1 + (Ch2-1)*nChannels);
                if(Ch1==Ch2)
                    plot(f,20*log10(abs(y(:,Ch1,Ch2))+eps),'r'); grid on; 
                     ylabel('psd (dB)'); xlabel('Frequency'); grid on; % power spectrum density
                else
                    confplot(f,y(:,Ch1,Ch2),yerr(:,Ch1,Ch2,1)); grid on; 
                     ylabel('chd (dB)'); grid on; % coherence density
                    
                     subplot(nChannels, nChannels, Ch2 + (Ch1-1)*nChannels);
                     
                    plot(f,unwrap(phi(:,Ch1,Ch2)),'k'); grid on; 
                    confplot(f,unwrap(phi(:,Ch1,Ch2)),unwrap(phierr(:,Ch1,Ch2)),unwrap(phierr(:,Ch1,Ch2)),'k'); grid on;
                    % unwrap means that the phase will be continuously
                    % reprensented: it avoids jumps; do mod for cyclic phase
                     ylabel('ph shift (rd)'); grid on; % phase shift
                end
            end
            ForAllSubplots('set(gca, ''FontSize'', 4)')
        end

        
plot(f,20*log10(abs(y(:,1,1))+eps),'r'); grid on; ylabel('psd (dB)'); xlabel('Frequency');  % power spectrum density



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% SPECTRUM COHERENCE  %%%%%%%%%%%%     %%%%%%%%%%%% MTCHGLONG %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


      %%%%%%%%%%%% compute %%%%%%%%%%%%

[y, f, t, phi] = mtchglong(EEG, nFFT, SR, WinLengthSample,[],3,'linear',[], FreqRange);

        %%%%%%%%%%%% plot %%%%%%%%%%%%  

 newplot; 
        for Ch1=1:nChannels, for Ch2 = Ch1:nChannels
                subplot(nChannels, nChannels, Ch1 + (Ch2-1)*nChannels);
                if Ch1==Ch2
                    if length(t)==1
                        
                        imagesc([0 1/f(2)],f,20*log10(abs(y(:,:,Ch1,Ch2))+eps)');
                        axis xy; colormap(jet); xlabel('Time'); ylabel('Frequency')
                    else
                        imagesc(t,f,20*log10(abs(y(:,:,Ch1,Ch2))+eps)');   % spectrogram
                        axis xy; colormap(jet); axis on;  colorbar; xlabel('Time'); ylabel('Frequency')
                    end

                else
                 
                    imagesc(t,f,(abs(y(:,:,Ch1,Ch2)))');    % coherogram
                    axis xy; colormap(jet);axis on; colorbar; xlabel('Time'); ylabel('Frequency')


                    subplot(nChannels, nChannels, Ch2 + (Ch1-1)*nChannels);
                    imagesc(t,f,squeeze(phi(:,:,Ch1,Ch2)+1.5)');    % phasogram
                    axis xy; colormap(jet);axis on; colorbar;xlabel('Time'); ylabel('Frequency')

                end
            end; end;
