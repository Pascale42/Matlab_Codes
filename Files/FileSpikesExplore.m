% stats = SpikeTrainExplore(Res, SampleRate, Epochs, Electrode, Options)
%
% plots various graphs pertaining to a spike train
%
% Epochs = time chunks within which spikes might have
% been collected (used for restricting to REM etc, also
% for calculating firing rates.
%
% if arg 1 is {Res,Clu}, it will loop thru all cells
% if it's {Res,Clu,CellNos} it will loop thru specified cells
% if it's {FileBase} it will load that.
% if it's {FileBase, CluNos} .. u get the idea
%
% also returns an array of statistics.
%
% Options can be 'nopause' in which case ... it won't pause.

function stats = SpikeTrainExplore(FileBase, varargin)
Par = LoadPar([FileBase '.par']);
SampleRate = 1e6./Par.SampleTime;
[fs, Epochs,Electrode, Options] = DefaultArgs(varargin, {SampleRate, [], [1:Par.nElecGps], 'nopause'});
GoodElectrodes=[];
for el=Electrode
	Clu = LoadClu([FileBase '.clu.' num2str(el)]);
	if max(Clu) > 1 
		GoodElectrodes = [GoodElectrodes, el];
	end
	
end	
[Res,Clu] = ReadEl4CCG(FileBase,GoodElectrodes);
%Res = load([FileBase '.res.' num2str(Electrode)]);
%Clu = LoadClu([FileBase '.clu.'  num2str(Electrode)]);
ToDo = 1:max(Clu);
for c=ToDo(:)'
    MyRes = Res(find(Clu==c));
    fprintf('Cell %d: firing rate %f\n', c, fs*length(MyRes)/(max(Res)-min(Res)));
    out = SpikeTrainExplore(Res(find(Clu==c)),fs,Epochs,Options);
    if ~isempty(out)
        stats(c) = out;
        fprintf('Press return.\n');
        if ~strcmp(Options, 'nopause')
            pause
        end
    end
end
return;
