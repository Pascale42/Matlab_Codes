function Trig2Evt(FileBase,varargin)

if FileExists('stlist')
	stfiles = LoadStringArray('stlist');
	GoThroughDb('Trig2Evt',stfiles,0,0,0);
else
	Par = LoadPar([FileBase '.par']);
	MakeEvtFile([FileBase '.trig'],[FileBase 'stim.evt'],'stim',1e6/Par.SampleTime);
end