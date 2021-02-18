function PlotSpectoSpd(FileBase, Channel, Window, FreqRange,fMode, Overwrite)



Par = LoadPar([FileBase '.xml']);
if isfield(Par,'lfpSampleRate')
    eSampleRate = Par.lfpSampleRate;
else
    eSampleRate = 1250;
end

switch fMode
    case 'compute'

        if ~isempty(Overwrite) | ~FileExists([FileBase '.spectro.mat'])
            % computes the spectrogram
            if FileExists([FileBase '.eeg'])
                try
                    Eeg = LoadBinary([FileBase '.eeg'], Channel,Par.nChannels)';
                catch
                    Eeg = LoadBinary([FileBase '.eeg'], Channel,Par.nChannels)';
                end
            else
                error('no eeg file! \n');
            end

            feeg = ButFilter(Eeg, 2, [FreqRange]/(eSampleRate/2),'bandpass');
            Spd = LoadEvt([FileBase '.spd.evt'],1250);

            %nFFT = 2^round(log2(2^11)); %compute nFFT according to different sampling rates
            SpecWindow = 2^round(log2(Window*eSampleRate));% choose window length as power of two
            nFFT = SpecWindow*4;
            weeg = WhitenSignal(Eeg,eSampleRate*2000,1);
            % [y,f,t]=mtcsglong(weeg,nFFT,eSampleRate,SpecWindow,[],2,'linear',[],FreqRange);
            save([FileBase '.spectro.mat'],'Eeg','feeg','Spd','Channels','FreqRange','Window');   % ,'-v6'
        else
            load([FileBase '.spectro.mat']);
        end

    case 'display'
        % now plots

        for i=1:length(Spd)
            subplot(311)
            x1=(Spd(i)-(Par.lfpSampleRate*1));
            x2=(Spd(i)+(Par.lfpSampleRate*5));
            [y,f,t]=mtcsglong(weeg(x1:x2),nFFT,eSampleRate,SpecWindow,[],2,'linear',[],FreqRange);
            imagesc(t,f,log(sq(y(:,:,1)))');axis xy; ylim([10 20]);
            subplot(312)
            plot(feeg(x1:x2))
            subplot(313)
            plot (Eeg(Spd(i)-(Par.lfpSampleRate*1):Spd(i)+(Par.lfpSampleRate*5)))
            waitforbuttonpress
        end

end



% for i=1:length(Spd)
%     subplot(311)
%     x1=((Spd(i)/Par.lfpSampleRate)-1);
%     x2=((Spd(i)/Par.lfpSampleRate)+5);
%     imagesc(t,f,log(sq(y(:,:,1)))');axis xy; ylim([10 20]); xlim([x1 x2]);
%     subplot(312)
%     plot(feeg(Spd(i)-(Par.lfpSampleRate*1):Spd(i)+(Par.lfpSampleRate*5)))
%     subplot(313)
%     plot (Eeg(Spd(i)-(Par.lfpSampleRate*1):Spd(i)+(Par.lfpSampleRate*5)))
%     waitforbuttonpress
% end



%
% subplot(211)
% imagesc(t,f,log(sq(y(:,:,1)))');axis xy; ylim([7 20]);
% subplot(212)
% plot (Eeg); axis xy
%
%
%
%  subplot(211)
% plot(feeg); axis xy
% subplot(212)
% plot (Eeg); axis xy

