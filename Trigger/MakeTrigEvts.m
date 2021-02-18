%function MakeTrigEvts(FileBase,Label)
function MakeTrigEvts(FileBase,varargin)

if FileExists('../stlabels')
	Labels = LoadStringArray('../stlabels');
	stlist = LoadStringArray('../stlist');
	ind = findcstr(stlist,FileBase,'strfind');
	DefLabel = Labels{ind};
else
	DefLabel = 'trig';
end
[Label]=DefaultArgs(varargin,{DefLabel});

Trig = load([FileBase '.trig']);

MakeEvtFile(Trig,[FileBase '.trig.evt'],1e6/Par.SampleTime,Label);


