function [ss, ds,fr] =PowerCoherence(filename)
%filename = '33705-6s';
par = LoadPar([filename '.par']);
eegpar = LoadEegPar(filename);

cxChan = find(strcmp(eegpar.ElecLoc,'c'));
cxChan = eegpar.ElecChannels{cxChan(1)}(1)+1;
hChan = find(strcmp(eegpar.ElecLoc,'h'));
hChan = eegpar.ElecChannels{hChan(1)}(1)+1;
cxeeg = readinsnglech([filename '.eeg'],par.nChannels,[cxChan]);
heeg = readinsnglech([filename '.eeg'],par.nChannels,[hChan]);
cxeeg = resample(cxeeg,1,5); 
cxsr = 250;
hsr=1250;
%rem = load([filename '.rem']);
%rem = round(rem/16);
%[eeg, ind] = SelectPeriods(eeg,rem,'c',0);
nTime = size(eeg,1);
deeg=[];seeg=[];spweeg=[];
win = sr*10;
nWins = ceil(nTime/win);
h=waitbar(0,'Wait..');
for i=1:nWins
    segwin = (1+(i-1)*win):min(i*win,nTime);
    %spindle and delta
    x = WaveletPowerRec(eeg(segwin,1),[1 4; 10 17],sr, 6);
    deeg = [deeg; x(:,1)];
    seeg = [seeg; x(:,2)];
    % spw power
    x = WaveletPowerRec(eeg(segwin,2),[150 220],sr, 6);
    spweeg = [spweeg;x];
   waitbar(i/nWins,h);
end
close(h);

h=waitbar(0,'Wait..');
for i=1:nWins
    segwin = (1+(i-1)*win):min(i*win,nTime);
    %spindle and delta
    x = WaveletPowerRec(eeg(segwin,1),[1 4; 10 17],sr, 6);
    deeg = [deeg; x(:,1)];
    seeg = [seeg; x(:,2)];
    % spw power
    x = WaveletPowerRec(eeg(segwin,2),[150 220],sr, 6);
    spweeg = [spweeg;x];
   waitbar(i/nWins,h);
end
close(h);


dfn = [filename '.dpow'];
bsave(dfn,dpow);
dfn = [filename '.spow'];
bsave(dfn,spow);
spwfn = [filename '.spwpow'];
bsave(spwfn,spwpow);
pow(:,1)=decimate(deeg,5);                      
pow(:,2)=decimate(seeg,1,5);
tmp=decimate(spweeg,1,5);
smwin =20;
tmp=conv(tmp,ones(smwin,1)/smwin);
pow(:,3)=tmp(ceil(smwin/2),end);
mtcsd(pow,2^14,125);
for i=1:9     
subplot(3,3,i)
xlim([0 1]);
end

% save(sfn,'spav','spstd','spw','f','t','s');
