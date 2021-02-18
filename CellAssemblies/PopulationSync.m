function PopulationSync(FileBase,varargin)

cxeEl = find(strcmp(eegpar.ElecLoc,'c'));
Par = LoadPar([FileBase '.par']);

[Res,Clu] = ReadEl4CCG(FileBase,cxEl);

States = {'REM','SWS'};

for i=1:2
	myState = States{i};
	myPer = SelectStates(FileBase,myState);
	
		[myRes 
SpikeCount(myRes, myClu, CluInd, BinSize,SampleRate,Norm)
