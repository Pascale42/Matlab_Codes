function FirRate = FileFiringRate(FileBase,varargin)
%function FirRate = FileFiringRate(FileBase,El)
El = DefaultArgs(varargin,{[]});

[Res,Clu,elclu] = ReadEl4CCG(FileBase,El);

Par = LoadPar([FileBase '.xml']);

FirRate = FiringRate(Res,Clu,[],Par.SampleRate);

