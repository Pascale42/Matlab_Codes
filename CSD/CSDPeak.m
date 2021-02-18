%function Peaks = CSDPeak(filename, Chans, FreqRange, Sign, Thr)
function Peaks = CSDPeak(filename,varargin)

[Chans, FreqRange, Sign, Thr] = DefaultArgs(varargin, {[3:5], [7 17], -1, 5});
csdChans = [Chans(1)-2:Chans(end)+2];
csd = MakeCSD(filename, csdChans, 2, 0);
csd = sum(csd,2);

%csd = FirFilter(csd, 20, FreqRange/625, 'bandpass');
csd = WaveletRec(csd,FreqRange,1250,6,1);
csd = NormMatrix(csd);
refr = 1250/(FreqRange(2));
Peaks = LocalMinimaK(-Sign*csd, refr,  -Thr);

figure
hist(csd,1000);
title(num2str(length(Peaks)));
if (nargout<0)
    msave([filename '.csdp' num2str(Chans(1)) '-' num2str(round(mean(FreqRange))) ], Peaks);
end

Peaks = [Peaks csd(Peaks)];