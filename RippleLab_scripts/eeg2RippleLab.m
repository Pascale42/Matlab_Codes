% function eeg2RippleLab
% 
% Loads given channels to RippleLab format
% 

function eeg2RippleLab

FileBase=gfb

nch=input('How many channels?');
channels=NaN(1,nch);
labels=cell(nch,1);
for n=1:nch
    labels(n,1)={input(['enter the name of the structure ' num2str(n) ': '], 's')};
    channels(n)= input(['enter the number for channel ' num2str(n) ': ']);
end

sts=input('Which state? : ', 's');
option=input('Option to add in the File name? : ', 's');

Data=GetEegState(FileBase, channels, 0, sts);
Data=Data./1000;

Header.Sampling=1250;
Header.Samples=repmat(size(Data,1),nch,1);
Header.Labels=labels;
Header.IniTime=[00,00,00];
Header.IniChannels=channels;

if strcmp(option, '')
    save([FileBase '.' sts '.mat'], 'Data', 'Header')
else
    save([FileBase '.' sts '.' option '.mat'], 'Data', 'Header')
end