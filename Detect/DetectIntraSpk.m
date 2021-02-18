function out = DetectIntraSpk(FileBase, varargin)
%function DetectIntraSpk(FileBase, SpkWidthMsec, CellRangeSec, IntraEl, Verbose, Refilter)
%e.g. FileBase ='8214-valentine.004.intra';
%     SpkWidthMsec = 5;
%     CellRangeSec = [4096 26406]; %seconds
%     IntraEl = 3 , where 3 is the electrode number for the
%     intra (if exists) and the channels for the current and voltage are
%     read from xml file, or [3 16 17], where 3 is the desired electrode
%     number for the intra to be, and 16 and 17 are respectively Vm and
%     Current channels counting from 0!!! 
Par = LoadPar([FileBase '.xml']);
SamplesInFile = FileLength([FileBase '.dat'])/2/Par.nChannels;

[SpkWidthMsec, CellRangeSec, IntraEl, Verbose, Refilter] = DefaultArgs(varargin,{[], [0 SamplesInFile/Par.SampleRate], Par.nElecGps, 0, 0});


CellRange = CellRangeSec*Par.SampleRate;
if CellRange(1)==0 CellRange(1)=1;end
BufferSize = Par.SampleRate*100;
Verbose = 0;
out.CellRangeSec = CellRangeSec;

%last spike group is intra. first channel - intra cell potential
if length(IntraEl)==1
    IntraChannel = Par.ElecGp{IntraEl}(1);
    CurrentChannel = Par.ElecGp{IntraEl}(2);
elseif length(IntraEl)==3
    tmp = IntraEl;
    IntraEl = tmp(1);
    IntraChannel = tmp(2);
    CurrentChannel = tmp(3);
else
    error('wrong IntraEl parameter');
end

out.IntraChannel = IntraChannel;
out.CurrentChannel = CurrentChannel;

if isempty(SpkWidthMsec)
    
   %estimate the spike width from the small chunk
   RefracPeriod = Par.SampleRate*2/1000; %assume 2 msec
   MedFilterLength = Par.SampleRate*15/1000; % assume 15 msec is max that can happen, be generous, in case of dendritic spike
   SegRange = CellRange(1)+[0 Par.SampleRate*10];
    dat_seg = LoadBinary([FileBase '.dat'], [IntraChannel CurrentChannel]+1, Par.nChannels, SegRange); 

    dcur_seg = abs(diff(dat_seg(:,2)));
    mdat_seg = dat_seg(:,1) - medfilt1(dat_seg(:,1),MedFilterLength, Par.SampleRate);
    lmc =  LocalMinima(-dcur_seg,100,std(dcur_seg));
    lmd =  LocalMinima(mdat_seg,RefracPeriod,0); %let's detect all of them and sort out later - then we don't need a threshold
    [lmd_art lmd_arti yi] = NearestNeighbour(lmd,lmc,'both', RefracPeriod);
    lmd(lmd_arti) = [];
    amp = mdat_seg(lmd);
    idx = kmeans(amp,2);
%     if Verbose
%         boxplot(amp,idx);
%     end
    idamp = accumarray(idx,amp,[2 1],@mean);
    [dummy goodid] = min(idamp);
    avspk = TriggeredAv(dat_seg(:,1), Par.SampleRate*2/1000, Par.SampleRate*5/1000, lmd(idx == goodid));
    [dummy spktr] = min(avspk);
    [dummy SpkWidth] = max(avspk(spktr:end));
    % if less than 32 or odd - correct
    SpkWidth = max(SpkWidth, round(Par.SampleRate*1.6/1000));
    SpkWidth = SpkWidth+mod(SpkWidth,2);
    out.AvSpk = avspk;
    out.t = [-Par.SampleRate*2/1000:Par.SampleRate*3/1000]/Par.SampleRate*1000;
    fprintf('estimated spike width is %d samples = %1.1f msec\n', SpkWidth, SpkWidth/Par.SampleRate*1000);
else
    SpkWidth = SpkWidthMsec*Par.SampleRate/1000;
end

MedFilterLength = SpkWidth*5;
Overlap = SpkWidth*10;
RefracPeriod = Par.SampleRate*1/1000; % 1 msec %SpkWidth;

out.SpkWidth =SpkWidth;
out.MedFilterLength = MedFilterLength;
out.RefracPeriod = RefracPeriod;

%extract intra channel and reverse to speedup median filter. 
if ~FileExists([FileBase '.idat'])
    fprintf('%s\n',cmd);
    cmd = sprintf('extractintra %s.dat %s.idat %d %d %d', FileBase, FileBase, Par.nChannels, IntraChannel, CurrentChannel);
    system(cmd);
end

%run median filter, which does: mdat = dat - medfilt1(dat,MedFilterLength, BufferSize);
%or  medianfilter -c 2 -s 2000000 -v -w 100 8214-valentine.004.intra.dat 8214-valentine.004.intra.mfil
if ~FileExists([FileBase '.mfil']) | Refilter
    fprintf('%s\n',cmd);
    cmd = sprintf('medianfilter -c %d -s %d -w %d %s.idat %s.mfil',2, BufferSize, MedFilterLength, FileBase, FileBase);
    system(cmd);
end

%%%%%%%%this section is to find the right threshold 
%cur =bload('8214-valentine.004.intra.dat',[1 inf], 2, 'short=>double',50*4+2);
%mdat = LoadBinary([FileBase '.mfil'],1,2,4,[],[],CellRange,100)';
cur = LoadBinary([FileBase '.dat'],CurrentChannel+1,Par.nChannels,CellRange,[], 100)';

dcur = abs(diff(cur));
[ch bh] = hist(dcur,10000);
lm = LocalMinima(ch,10);
[mch mchi] = max(ch);
dcurThresh = lm(lm>mchi);
dcurThresh = bh(dcurThresh(1));
% minimal current step is 0.1 nA, which assuming 10V range is ~ 2000
% so we require the threshold to be minimum 1000
dcurThresh = max(dcurThresh, 1000);
out.dcurThresh = dcurThresh;

%lm_dcur = LocalMinima(-dcur,100,-dcurThresh);

% mdat(mdat>0) = 0;
% absmdat = abs(mdat);
% %[dh bh1] = hist(absmdat,50);
% lm_absmdat = LocalMinima( -absmdat(absmdat<prctile(absmdat,95)), RefracPeriod, 0);
% [dh bh1] = hist(absmdat(lm_absmdat),100);
% 
% lm = LocalMinima(dh,10);
% [mch mchi] = max(dh);
% mdatThresh = lm(lm>mchi);
% mdatThresh = bh1(mdatThresh(1));
% lm_mdat = LocalMinima(mdat,10,-mdatThresh);

nT = diff(CellRange)+1;
Step = BufferSize-Overlap;
ns = ceil((nT-BufferSize)/(BufferSize-Overlap)+1);
Res = [];
Clu = [];
Amp = [];
for s=1:ns
    SegBeg =CellRange(1)+(s-1)*(BufferSize-Overlap)+1; 
    SegEnd =min(SegBeg+BufferSize-1,CellRange(1)+nT-1);
    
    SegRange = [SegBeg SegEnd];
    if s>1 & s<ns
        UseSegRange = [Overlap+1 BufferSize-Overlap];
    elseif s==1
        UseSegRange = [1 BufferSize-Overlap];
    elseif s==ns
        UseSegRange = [Overlap+1 diff(SegRange)+1];
    end
    
    mdat_seg = LoadBinary([FileBase '.mfil'], 1, 2, 4, [],[],SegRange)'; 
    cur_seg = LoadBinary([FileBase '.idat'], 2, 2, 4, [],[],SegRange)'; 
    dcur_seg = abs(diff(cur_seg));

    lmc =  LocalMinima(-dcur_seg,100,-dcurThresh);

    lmd =  LocalMinima(mdat_seg,RefracPeriod,0); %let's detect all of them and sort out later - then we don't need a threshold
    
    [lmd_art lmd_arti yi] = NearestNeighbour(lmd,lmc,'both', RefracPeriod);
    
    %lmd(lmd_arti) = [];
    lmd = SelectPeriods(lmd,UseSegRange,'d',1);
    myClu = 2*ones(size(lmd));
    myClu(lmd_arti) = 0;
    Clu = [Clu; myClu];
    Res = [Res; lmd+SegBeg];
    Amp = [Amp; mdat_seg(lmd)];
    if length(Amp)~=length(Clu) break;end
end

%now threshold the spikes by Amp 
GoodAmp = Amp(Clu==2);
idx = kmeans(GoodAmp,2);

mAmp = accumarray(idx,GoodAmp,[2 1],@mean);
if mAmp(1)>mAmp(2)
    tmp=idx; tmp(idx==1)=2; tmp(idx==2)=1; 
    [mAmp(1) mAmp(2)] = deal(mAmp(2), mAmp(1));
end
[AmpDens b] = ksdensity(GoodAmp,'npoints',1000);
[dummy mAmpInd]= NearestNeighbour(b, mAmp, 'both');
[dummy minDensInd] = min(AmpDens(mAmpInd(1):mAmpInd(2)));
densThresh = mAmp(1)+b(minDensInd);

% [mch mchi] = max(h);
% lm = LocalMinima(h,5);
% histThresh = lm(lm<mchi);
% histThresh = b(histThresh(end));

AmpThr(1) = densThresh;
AmpThr(2) = max(GoodAmp(idx==goodid));
AmpThr(3) = mean(mAmp);
%now do the same in the running window
%later!!!!

AmpThr = max(AmpThr); %now we chose the threshold 
out.AmpThr = AmpThr;

%Res = Res(Amp<AmpThr);
Clu(Amp>AmpThr)=1;

nRes = length(Res);

%before we write all new files, let's backup old ones

cmd = sprintf('mv %s.res.%d %s.res.%d.old',FileBase, IntraEl, FileBase, IntraEl);
fprintf('%s\n',cmd);
system(cmd);
cmd = sprintf('mv %s.clu.%d %s.clu.%d.old',FileBase, IntraEl, FileBase, IntraEl);
fprintf('%s\n',cmd);system(cmd);
cmd = sprintf('mv %s.spk.%d %s.spk.%d.old',FileBase, IntraEl, FileBase, IntraEl);
fprintf('%s\n',cmd);system(cmd);
cmd = sprintf('mv %s.m1m2.%d %s.m1m2.%d.old',FileBase, IntraEl, FileBase, IntraEl);
fprintf('%s\n',cmd);system(cmd);
cmd = sprintf('mv %s.fet.%d %s.fet.%d.old',FileBase, IntraEl, FileBase, IntraEl);
fprintf('%s\n',cmd);system(cmd);
cmd = sprintf('mv %s.mm.%d %s.mm.%d.old',FileBase, IntraEl, FileBase, IntraEl);
fprintf('%s\n',cmd);system(cmd);
cmd = sprintf('cp %s.xml %s.xml.old',FileBase, FileBase);
fprintf('%s\n',cmd);system(cmd);



msave([FileBase '.res.' num2str(IntraEl)],Res);
msave([FileBase '.clu.' num2str(IntraEl)],Clu);

%process_getseg $base.res.$electrodeGroup $base.fil $base.spk.$electrodeGroup -n $nChannels -c $electrodeList -s $samplesInWaveform -p $peakSample
cmd = sprintf('process_getseg %s.res.%d %s.mfil %s.spk.%d -n 2 -c 0,1 -s %d -p %d',...
    FileBase, IntraEl, FileBase, FileBase, IntraEl, SpkWidth, round(SpkWidth/3));
system(cmd);

%process_variance $base.fil $base.res.$2 $base.m1m2.$electrodeGroup -o $offset -n $nChannels -c $electrodeList -s $samplesInWaveform -p $peakSample
cmd = sprintf('process_variance %s.mfil %s.res.%d  %s.m1m2.%d -o 0 -n 2 -c 0,1 -s %d -p %d',...
    FileBase, FileBase, IntraEl, FileBase, IntraEl, SpkWidth, round(SpkWidth/3));
system(cmd);

beforePeak = round(SpkWidth/3)-5;
afterPeak = round(SpkWidth/3)-5;
nSamplesInPCA = beforePeak + afterPeak;

%process_feature $base.spk.$electrodeGroup $base.m1m2.$electrodeGroup $base.fet.$electrodeGroup $base.mm.$electrodeGroup -o $offset -n $nElectrodesInGroup
% -s $samplesInWaveform -p $peakSample -b $beforePeak -a $afterPeak -t "eigen" -f $nFeatures -v $nSamplesInPCA
cmd = sprintf('process_feature %s.spk.%d %s.m1m2.%d  %s.fet.%d %s.mm.%d -o 0 -n %d -s %d -p %d -b %d -a %d -t "eigen" -f 5 -v %d',...
    FileBase, IntraEl, FileBase, IntraEl, FileBase, IntraEl, FileBase, IntraEl, 2, SpkWidth,round(SpkWidth/3), beforePeak, afterPeak, nSamplesInPCA);
fprintf([cmd '\n']);
system(cmd);
 
 
%#Add the content of the res file (spike time) as the last column of the fet file
%#Add the number of features in the first line of the feature file
%head -1 $base.fet.$electrodeGroup | awk '{print $1+1}'
%>jnk.$electrodeGroup
cmd = sprintf('head -1 %s.fet.%d | awk ''{print $1+1}'' >jnk', FileBase, IntraEl);
system(cmd);

%tail +2 $base.fet.$electrodeGroup | paste -d " " -
%$base.res.$electrodeGroup >> jnk.$electrodeGroup
cmd = sprintf('tail +2 %s.fet.%d | paste -d " " - %s.res.%d >> jnk', FileBase, IntraEl,FileBase, IntraEl); 
system(cmd);

%mv jnk.$electrodeGroup $base.fet.$electrodeGroup
cmd = sprintf('mv jnk %s.fet.%d', FileBase, IntraEl);
system(cmd);

%#Keep only the first line ?
%tail -1 $base.mm.$electrodeGroup >> $base.mm.$electrodeGroup
cmd = sprintf('tail -1 %s.mm.%d >> %s.mm.%d', FileBase, IntraEl,FileBase, IntraEl);
system(cmd);

system('rm -f  jnk');

%update the fileds in spike group in xml file

%quick hach assuming that format of xml file is not changed.

xml = xmltools([FileBase '.xml']);
%find the chgrp tag
for l=1:length(xml.child(2).child)
    if strcmp(xml.child(2).child(l).tag, 'spikeDetection')
        chgrpid = l;
    end
end
chgrps =  xml.child(2).child(chgrpid).child(1);
mygrp = chgrps.child(IntraEl);
mygrp.child(1).child(1) = struct('tag','channel','attribs',struct('name','','value',''),'value',num2str(IntraChannel),'child',[]);
mygrp.child(1).child(2) = struct('tag','channel','attribs',struct('name','','value',''),'value',num2str(CurrentChannel),'child',[]);
mygrp.child(2) = struct('tag','nSamples','attribs',struct('name','','value',''),'value',num2str(SpkWidth),'child',[]);
mygrp.child(3) = struct('tag','peakSampleIndex','attribs',struct('name','','value',''),'value',num2str(round(SpkWidth/3)),'child',[]);
mygrp.child(4) = struct('tag','nFeatures','attribs',struct('name','','value',''),'value','5','child',[]);

xml.child(2).child(chgrpid).child(1).child(IntraEl) = mygrp;

xmltools(xml,[FileBase '.xml']);

