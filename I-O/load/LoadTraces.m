function LoadTraces(FileBase,Channels,

Par = LoadPar([FileBase '.xml']);

eeg  = LoadBinary([FileBase '.eeg'], Channels, Par.nChannels);
vm = MakeVm(FileBase,'get');