%function Rips = DetectRipples(FileBase, varargin:: Channels, FreqRange, State, Threshold,  Mode, Overwrite)
% gives out the structure: Rip.t - time of the ripple (eeg sampling rate),
% Rip.len = length is msec, Rip.pow - power of the ripple


function Rips = DetectRipples(FileBase, varargin)

[Channels, FreqRange, State, Threshold,  Mode, Overwrite] = DefaultArgs(varargin,{[],[100 250], [], 5, 'long',1});
if isempty(Channels)
%     if FileExists([FileBase '.eegseg.par'])
%         Channels = load([FileBase '.eegseg.par']);
%         Channels = Channels(1)+1;
%     else
        error('Channels not specified');
%     end
end

if ~FileExists([FileBase '.spw']) || Overwrite
    
    Par = LoadPar([FileBase '.xml']);
    eSampleRate = Par.lfpSampleRate;
    spw = sdetect_a([FileBase '.eeg'],Par.nChannels,Channels,FreqRange(end),FreqRange(1),Threshold, 0, eSampleRate);

    if strcmp(Mode,'long') && ~isempty(spw)

        seglen = max(spw(:,3)-spw(:,1));
        seglen = seglen + mod(seglen,2);
        eeg = LoadBinary([FileBase '.eeg'],Channels(1),Par.nChannels);
        [seg, complete] = GetSegs(eeg,spw(:,2)-seglen/2,seglen,[]);
        seg = squeeze(seg);
        pow = FirFilter(seg,2,[120 230]/(eSampleRate/2),'bandpass');
        Rips.pow = mean(abs(pow),1)';

        Rips.t = spw(complete,2);
        Rips.len = (spw(complete,3)-spw(complete,1))*1000./eSampleRate;
        
        % ajout du Avril15th 2015
    if strcmp(State, 'SWS')
        SWS=load([FileBase '.sts.SWS']);
    [Rips.t, ind]=SelPerDiscr(Rips.t,SWS,1,1);
%         [Rips.T, ind]=SelPerDiscr2(Rips.t,SWS);
    Rips.len=Rips.len(ind);
    Rips.pow=Rips.pow(ind);
    MakeEvtFile(Rips.t(:), [FileBase '.spw.evt'],'spw',eSampleRate,1);
    end
        
%         msave([FileBase '.spw'],[Rips.t Rips.pow(:) Rips.len(:)]);   changed April15th 2015
        save([FileBase '.spw.mat'],'Rips', 'Threshold');
         MakeEvtFile(Rips.t(:), [FileBase '.rip.evt'],'rip',eSampleRate,1);

        
    end

end
