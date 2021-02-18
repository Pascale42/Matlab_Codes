
rep=1;
if (rep==0)
    filebase='836602s';
    
    Par = LoadPar([filebase '.par']);
    EegPar = LoadEegPar(filebase);
    
    %getting just first cortex or hippocampal electrode ..diregarding if it is best or not..
    cxeleeg=findcell(EegPar.ElecLoc,'c','eq'); cxeleeg=cxeleeg(1); % for eeg selection
    hipeleeg=findcell(EegPar.ElecLoc,'h','eq'); hipeleeg=hipeleeg(1);
    
    %for units
    cxel=findcell(Par.ElecLoc,'c','eq'); cxel=cxel(1);
    hipel=findcell(Par.ElecLoc,'h','eq'); hipel=hipel(1);
    
    %reading cx eeg form 1st channel in cxel
    eegcx=readsinglech([filebase '.eeg'], Par.nChannels, EegPar.ElecChannels{cxeleeg}(1)+1);
    eeghip=readsinglech([filebase '.eeg'], Par.nChannels, EegPar.ElecChannels{hipeleeg}(1)+1);
    
    eegcxdelt = filtereeg(eegcx, [1.5 4.5], [0.5 5.5], 1250);
    cxres = load([filebase '.res.' num2str(cxel)]);
    hipres = load([filebase '.res.' num2str(hipel)]);
end
hipres = load([filebase '.res.' num2str(hipel)]);
resamplerate = 50;
eegcx1 = resample(eegcx, 1, resamplerate);

cxresb = 1 + floor( (cxres-1)/16/resamplerate ); 
cxreseeg = Accumulate(cxresb, 1, length(eegcx1));
cxreseeg0 = cxreseeg - mean(cxreseeg);


% figure(1)
% Comodugram([eegcx1, cxreseeg], 2^8, 1250/resamplerate, [0 100]);
% 
% figure(2)
% Comodugram2([eegcx1, cxreseeg], 2^8, 1250/resamplerate, [0 100]);

hipresb = 1 + floor( (hipres-1)/16/resamplerate); 
hipreseeg = Accumulate(hipresb, 1, length(eegcx1));
hipreseeg0 = hipreseeg - mean(hipreseeg);


figure(1)
Comodugram([eegcx1, hipreseeg], 2^8, 1250/resamplerate, [0 100]);


