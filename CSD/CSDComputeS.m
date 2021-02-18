%function PowerComputeS(filename,ch)
function PowerComputeS(filename,varargin)
ch = DefaultArgs(varargin,{[5 9]});
%filename = '33705-6s';
par = LoadPar([filename '.par']);
% if (nargin<2)
%     eegpar = LoadEegPar(filename);
% 
% 	cxChan = find(strcmp(eegpar.ElecLoc,'c'));
% 	cxChan = eegpar.ElecChannels{cxChan(1)}(1)+1;
% 	hChan = find(strcmp(eegpar.ElecLoc,'h'));
% 	hChan = eegpar.ElecChannels{hChan(1)}(1)+1;
% else
    cxChans=ch(1)+[-2 0 2];
    cxChand=ch(2)+[-2 0 2] ;

% end
eeg = readmulti([filename '.eeg'],par.nChannels,[cxChans cxChand ]);
%eeg = resample(eeg,1,5); 

csd = eeg(
