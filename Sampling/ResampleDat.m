function ResampleDat(FileBase,varargin)

Par = LoadPar([FileBase '.par']);
tfs =20000;
fs = 1e6/Par.SampleTime;

%greatest common divisor
cdiv = gcd(tfs,fs);

resamplup = tfs/cdiv;

resampldown = fs/cdiv;

ResampleFile([FileBase '.dat'],[FileBase '.cnv'],Par.nChannels,resampldown, resamplup);

