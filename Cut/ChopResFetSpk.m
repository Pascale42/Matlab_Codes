%
% function ChopResFetSpk(FileBase, eldgrp, Evt, cuttime);
%
% Chop res fet and spk according to time in evt file
%
% Evt='cut'; chop the file according to begining & end of FileBase.evt.cut file
% created in neuroscope; tag must be called 'cut'
% Evt='cel'; chop the file according to begining & end of FileBase.evt.cel file
%   cuttime=1 if wants to chop after 1 time point (begining of evt.cel file)
%   cuttime=2 if wants to chop in between 2 points (begining & end of evt.cel file)

function ChopResFetSpk(FileBase, eldgrp, Evt, cuttime);

if ~FileExists([FileBase '.res.' eldgrp '.big']);

    % Rename fet, res and spk
    s = sprintf('mv %s.fet.%d %s.fet%d.big', FileBase, eldgrp, FileBase, eldgrp);
    fprintf('%s\n',s);
    system(s);

    s = sprintf('mv %s.res.%d %s.res.%d.big', FileBase, eldgrp, FileBase, eldgrp);
    fprintf('%s\n',s);
    system(s);

    s = sprintf('mv %s.spk.%d %s.spk.%d.big', FileBase, eldgrp, FileBase, eldgrp);
    fprintf('%s\n',s);
    system(s);
else
end


% Loading
par=LoadPar([FileBase '.xml']);

fid = fopen([FileBase '.res.' eldgrp '.big'], 'r');
t = fscanf(fid,'%d');
fclose(fid);

switch Evt
    case 'cut'
        [EvtRes] = LoadEvt([FileBase '.evt.cut'],par.SampleRate,'cut');
        r = find(t>=EvtRes(1) & t<=EvtRes(end));

    case 'cel'

        [EvtRes] = LoadEvt([FileBase '.evt.cel'],par.SampleRate,'cel');

        if cuttime=1
            r = find(t>=EvtRes(1));
        else
            r = find(t>=EvtRes(1) & t<=EvtRes(end));
        end
end

T = t(r);

% creates a new res file
msave([FileBase '.res.' eldgrp], T);

% chop and makes fet file
TailFet = length(t)-r(1)+1;
HeadFet = r(end)-r(1)+1;

sfet = sprintf('head -n 1 %s.fet.%d.big > %s.fet.%d', FileBase, eldgrp, FileBase, eldgrp);
fprintf('%s\n',sfet);
system(sfet);
sfet = sprintf('tail -n %d %s.fet.%d.big | head -n %d >> %s.fet.%d', TailFet, FileBase, eldgrp, HeadFet, FileBase, eldgrp);
fprintf('%s\n',sfet);
system(sfet);

% chop and makes spk file
    % # bytes for one spike (=one block)
    % By = par.SpkGrps(1,5).nSamples * 2 * par.nChannels; not par.nChannels!
    % nChannels in electrode group!
By = par.SpkGrps(1,eldgrp).nSamples * 2 * length(par.SpkGrps(1,eldgrp).Channels);
TailSpk = TailFet * By;
HeadSpk = HeadFet * By;

sspk = sprintf('tail -c %d %s.spk.%d.big | head -c %d > %s.spk.%d', TailSpk, FileBase, eldgrp, HeadSpk, FileBase, eldgrp);
fprintf('%s\n',sspk);
system(sspk);
