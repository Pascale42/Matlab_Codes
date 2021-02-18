%Triggers = function DetectTrigger(FileBase, TrigChan, Sign, FileToUse, Control, nChannels, SampleRate)
% detects the trigger signal from the TrigChan (by default last)
% on the downstroke of the negative step according to Sign (default -1)
% Trigger in 'eeg' sampling units 
function Triggers = DetectTrigger(FileBase, varargin)

Par= LoadPar([FileBase '.par']);
[TrigChan, Sign  , FileToUse, Control, nChannels,SampleRate] = DefaultArgs(varargin, {Par.nChannels, -1, 'eeg',0,Par.nChannels,[] });
if strcmp(FileToUse,'eeg') & isempty(SampleRate) 
    SampleRate = 1250;
elseif isempty(SampleRate) 
    SampleRate = 1e6/Par.SampleTime;
end
%SampleRate
dat = readsinglech([FileBase '.' FileToUse], nChannels, TrigChan);

dat = dat*Sign;
Thr = mean(dat)+10*std(dat);
Triggers = SchmittTrigger(dat,Thr, Thr);


nTrigs = length(Triggers);
TrigInterval = diff(Triggers)/SampleRate;
[TrigWave, Complete] = GetSegs(dat, Triggers, round(1.250*100));

[u,s,v] = svd(TrigWave(:,Complete),0);
pca1 = abs(v(:,1));


figure(445)
subplot(311)
plot(unity(dat))
hold on
Lines(Triggers, [],' k');

subplot(312)
plot(TrigWave)
title('Trigger waves');
subplot(325)
hist(pca1);
title('PCA1 of the trigger step');
subplot(326)
hist(TrigInterval);
title('Inter trigger interval');
while Control
    Control = input('Do you want to change (0/1)?');
    if Control
        figure(445)
        subplot(325)
        [PcaThr, dummy ] = ginput(1);
        GoodTrig = find(pca1>PcaThr);
        BadTrig = setdiff([1:nTrigs],GoodTrig);
        subplot(311)
        Lines(Triggers(GoodTrig), [],' r');
        
        subplot(312)
        cla
        plot(TrigWave(:,GoodTrig),'r');
        hold on
        plot(TrigWave(:,BadTrig),'k');
                
        subplot(325)
        cla
        bins = linspace(min(pca1)-eps, max(pca1)+eps,10);
        hb = hist(abs(v(BadTrig,1)),bins);
        
        hg = hist(abs(v(GoodTrig,1)),bins);
        bar(bins,hg,'r');
        hold on
        bar(bins,hb,'k');
        legend({'good','bad'});
        title('PCA1 of the trigger step');

        subplot(326)
        bins = linspace(min(TrigInterval),max(TrigInterval),10);
        hg = hist(diff(Triggers(GoodTrig))/SampleRate, bins);
        hb = hist(diff(Triggers(BadTrig))/SampleRate, bins);
        bar(bins,hg,'r');
        hold on
        bar(bins,hb,'k');
        legend({'good','bad'});
        title('Inter trigger interval');
        
    end       
end

fprintf('Detected %d triggers\n',length(Triggers));
pause
close(445)
if SampleRate ==1250
    Triggers = Triggers * 16;
end
msave([FileBase '.trig'], Triggers);
fprintf('Triggers in raw data sampling rate are saved into file %s.trig', FileBase);