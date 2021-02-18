%function histogramme1 : 
%FileName
%anatgrps
%channels
%Threshold1
%Threshold2
%anatgrpstest
%channelstest
%channelstestfin
%maxmean=1 -> max

function out=histogramme1(FileName,anatgrps,channels,varargin) %add anatgrpstest, channelstest

[Threshold1,Threshold2,anatgrpstest, channelstest, channelstestfin, maxmean]...
    = DefaultArgs(varargin,{5,2,[],[],[],[]});

%maxmean=1->max
%maxmean=0->mean
    
if anatgrps==anatgrpstest & channels==channelstest
    error('Channel test = Channel');
end 

if isempty(channelstestfin)       %si channelstestfin est vide
    if isempty(anatgrpstest) && isempty(channelstest)    %si anatgrpstest et channelstest vides
    FileBase=FileName;
    Par= LoadPar([FileBase '.xml']);
    Channels=Par.AnatGrps(anatgrps).Channels([channels])+1;
    %Channelstest=Par.AnatGrps(anatgrpstest).Channels([channelstest])+1;
    %EEG=LoadBinary([FileBase '.eeg'],Channels, Par.nChannels);
    out.Ripples=DetectRipples1(FileBase, Channels,[],[100 200],[Threshold1 Threshold2]);

    
    else
    FileBase=FileName;
    Par=LoadPar([FileBase '.xml']);
    Channels=Par.AnatGrps(anatgrps).Channels([channels])+1;
    Channelstest=Par.AnatGrps(anatgrpstest).Channels([channelstest])+1;
    %EEG=LoadBinary([FileBase '.eeg'],Channels, Par.nChannels);
    out.Ripples=DetectRipples1(FileBase, Channels,Channelstest,[100 200],[Threshold1 Threshold2]);
   
    end
else
    for i=channelstest:1:channelstestfin

    FileBase=FileName;
    Par=LoadPar([FileBase '.xml']);
    Channels=Par.AnatGrps(anatgrps).Channels([channels])+1;
    Channelstest=Par.AnatGrps(anatgrpstest).Channels([i])+1;
    %EEG=LoadBinary([FileBase '.eeg'],Channels, Par.nChannels);
    out.Ripples=DetectRipples1(FileBase,Channels,Channelstest,[100 200],[Threshold1 Threshold2]);
    
    out.nb(i-channelstest+1)=length(out.Ripples.t);
    
     if 1
      bar(out.nb);
     end 
    end
if maxmean 
        [m,j]=max(out.nb);
        Channelstest=Par.AnatGrps(anatgrpstest).Channels([j])+1;
        out.Ripples=DetectRipples1(FileBase,Channels,Channelstest,[100 200],[Threshold1 Threshold2]);
        
    else
        moy=mean(out.nb);
        nb1=abs(out.nb-moy);
        [m, j]=min(nb1);
        Channelstest=Par.AnatGrps(anatgrpstest).Channels([j])+1;
        out.Ripples=DetectRipples1(FileBase,Channels,Channelstest,[100 200],[Threshold1 Threshold2]);
       
        
    end
end
    out.FileBase=FileBase;
    out.Par=Par; 
    %out.Channels=Channels;
    %out.Channelstest=Channelstest;
    %out.EEG=EEG;
    out.Threshold=[Threshold1 Threshold2];
end 

function Rips = DetectRipples1(FileBase, varargin)
[Channels, Channelstest, FreqRange, Threshold,  Overwrite] = DefaultArgs(varargin,{[],[],[100 250],[5 2],1}); %addition Channelstest
if isempty(Channels)
    if FileExists([FileBase '.eegseg.par'])
        Channels = load([FileBase '.eegseg.par']);
        Channels = Channels(1)+1;
    else
        error('Channels not specified'); 
    end
end

% addition
if isempty(Channelstest) 
   
    Par = LoadPar([FileBase '.xml']);
    eSampleRate = Par.lfpSampleRate;
    
    lfp = LoadBinary([FileBase '.eeg'],Channels(1),Par.nChannels,4)'; 
    
    Rips = DetectOscillations(lfp, FreqRange, 5, eSampleRate, Threshold);
    %TODO: get Oscillation Frequency for each ripple
    msave([FileBase '.spw'],[Rips.t Rips.pow(:) Rips.len(:)]);
    MakeEvtFile(Rips.t(:), [FileBase '.rip.evt'],'rip',eSampleRate,1);
    save([FileBase '.' mfilename '.mat'],'Rips'); 
    
%addition
else if ~FileExists([FileBase '.spw']) | Overwrite
        
    Par = LoadPar([FileBase '.xml']);
    eSampleRate = Par.lfpSampleRate;
    
    lfp = LoadBinary([FileBase '.eeg'],Channels(1),Par.nChannels,4)'; 
    lfpt = LoadBinary([FileBase '.eeg'],Channelstest(1),Par.nChannels,4)'; %addition lfpt
    
    lfps=lfp-lfpt;
    
    Rips = DetectOscillations(lfps, FreqRange, 5, eSampleRate, Threshold);
    %TODO: get Oscillation Frequency for each ripple
    msave([FileBase '.spw'],[Rips.t Rips.pow(:) Rips.len(:)]);
    MakeEvtFile(Rips.t(:), [FileBase '.rip.evt'],'rip',eSampleRate,1);
    save([FileBase '.' mfilename '.mat'],'Rips');
    
    end 
end

end

function sdetected = sdetect_a(inname,numchannel,channels,highband,lowband,thresholdf,verbose,sampl)

% low threshold version

%%% parameters to play with

% parameters for program flow control

showfiltchar = 0;
showrealanal = 0;
%verbose = 1;
pause off;

% parameters for detection tuning

%sampl = 1250;      % sampling rate of the eeg file in Hz
if nargin<4
    highband = 230;    % bandss filter range (230Hz to 100Hz)
    lowband = 140;     %
    thresholdf=10.0;
end
forder = 100;  % filter order has to be even; .. the longer the more
% selective, but the operation will be linearly slower
% to the filter order

avgfilorder = 21; % do not change this... length of averaging
% filter

avgfilterdelay = floor(avgfilorder/2);  % compensated delay period

forder = ceil(forder/2)*2;           %make sure filter order is even

buffersize = 2^24; % buffer size ... the larger the faster, as long as
% (real) memory can cover it.


min_sw_period = 30 ; % minimum sharpwave period = 25ms ~ 6 cycles
% of ripples

max_sw_period = 150; % maximum sharpwave period = 150ms ~ 30 cycles
% of ripples

min_isw_period = 5; % minimum inter-sharpwave period = 30 ms;


%thresholdf = 10.0;   % threshold for ripple detection
% based upon the standard deviation measure
if isstr(inname)
    datafile = fopen(inname,'r'); % specify input and output file
end
regionexpand = 10;  % add this number of samples in the end of
% detection and subtract this number of samples
% in the beginning of the detection. 10->8ms.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% changing something below this line will result in changing the
% algorithm (i.e. not the parameter)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% %
% calculate the convolution function (passbands are normalized to
% the Nyquist frequency

firfiltb =  fir1(forder,[lowband/sampl*2,highband/sampl*2]);
avgfiltb = ones(avgfilorder,1)/avgfilorder;

% determine the threshold
sdat = [];
%verdisp('computing threshold...',verbose);
thresholdbuffersize = 2^14;
while (~feof(datafile))
    [thresholdbuffer, count] = fread(datafile,[numchannel,thresholdbuffersize],'int16');
    thresholdbuffer=thresholdbuffer';
    if ~isempty(thresholdbuffer)
        initvec = thresholdbuffer(1,:);
        [thresholdbuffer,initvec] = adjovfl(thresholdbuffer,initvec);
        filtered_data = filter(firfiltb,1,thresholdbuffer(:,channels));
        filtered_data = abs(filtered_data);
        filtered_data2 = filter(avgfiltb,1,filtered_data);
        sdat = [sdat;sum(filtered_data2,2)];
    end
end

[thresholdbuffer, count] = fread(datafile,[numchannel,thresholdbuffersize],'int16');

sthreshold = std(sdat)*thresholdf;
clear thresholdbuffer;
clear sdat;
frewind(datafile);
%verdisp('done',verbose);
pause;

% overlap buffer length;
overlap = forder;
overlap1 = overlap+1;
overlap11 = overlap-1;

overlap2 = overlap/2;
overlap21 = overlap2+1;

filcond = [];
sdetected = [];
offset = 0;
threshstate = 0;

% initial overlap buffer ... which is actually the mirror image
% of the initial portion ... make sure that you 'rewind' the file
%verdisp('detecting ...',verbose);

overlapbuffer = fread(datafile,[numchannel,overlap/2],'int16');
frewind(datafile);
overlapbuffer = transpose(fliplr(overlapbuffer));
[datasegment,count] = fread(datafile,[numchannel,buffersize],'int16');
datasegment2 = [overlapbuffer;datasegment'];

initvec = datasegment2(1,:);
[datasegment2,initvec] = adjovfl(datasegment2,initvec);

filtered_data = filter(firfiltb,1,datasegment2(:,channels));
filtered_data = abs(filtered_data);
[filtered_data2 ,filcond] = filter(avgfiltb,1,filtered_data, ...
    filcond);
filtered_data3 = sum(filtered_data2,2);

[dd,threshstate] = threshres(filtered_data3(overlap1:end), ...
    sthreshold,threshstate);
sdetected = dd;

offset = offset + length(filtered_data(overlap1:end,1));

% show the frequency and phase response
[h,w] =  freqz(firfiltb,1,buffersize,sampl);

if showfiltchar
    figure(1);
    clf;
    w1 = subplot(2,2,1);plot(firfiltb);
    title('filter function');
    xlabel('time');
    set(w1,'XLim',[1,forder]);
    w2 = subplot(2,2,2);plot(w,unwrap(angle(h)));
    title('phase');
    xlabel('frequency(Hz)');
    ylabel('radian');
    set(w2,'XLim',[0,sampl/2]);
    w3 = subplot(2,2,3);plot(w,abs(h));
    title('amplitude');
    xlabel('frequency(Hz)');
    set(w3,'XLim',[0,sampl/2]);
    w4 = subplot(2,2,4);plot(w,log10(abs(h))*20);
    title('amplitude in decibels');
    xlabel('frequency(Hz)');
    set(w4,'XLim',[0,sampl/2]);
end

if showrealanal
    figure(2)
    plotbuffer = buffersize/2;
    w5 = subplot(2,1,1);
    plot([1:plotbuffer]/sampl,datasegment(1,1:plotbuffer));
    title('sample raw data');
    xlabel('time(s)');
    axis([0,plotbuffer/sampl,min(datasegment(1,1:plotbuffer)), ...
        max(datasegment(1,1:plotbuffer))]);
    w6 = subplot(2,1,2);
    plot([1:plotbuffer]/sampl, ...
        filtered_data(overlap1:plotbuffer+overlap,1));
    title('sample filtered data');
    xlabel('time(s)');
    axis([0,plotbuffer/sampl,...
        min(filtered_data(overlap1:plotbuffer+overlap,1)),...
        max(filtered_data(overlap1:plotbuffer+overlap,1))]);
end

% do the rest (this routine can be made faster, if you somehow

while ~feof(datafile),
    fseek(datafile,-2*numchannel*overlap,0);
    datasegment = fread(datafile,[numchannel,buffersize],'int16')';
    
    [datasegment,initvec] = adjovfl(datasegment,initvec);
    
    filtered_data = filter(firfiltb,1,datasegment(:,channels));
    filtered_data = abs(filtered_data);
    [filtered_data2 ,filcond] = filter(avgfiltb,1,...
        filtered_data,filcond);
    filtered_data3 = sum(filtered_data2,2);
    [dd,threshstate] = threshres(filtered_data3(overlap1:end,1), ...
        sthreshold,threshstate);
    
    dd = (abs(dd) + offset) .* sign(dd);
    sdetected = [sdetected;dd];
    
    if showrealanal
        axes(w5);
        plot(datasegment(overlap2:end-overlap2,1));
        axes(w6);
        hold off
        plot(filtered_data3(overlap1:end),'b');
        hold on;
        plot([1,buffersize],[sthreshold,sthreshold],'r:');
        plot(abs(dd)-offset,ones(length(dd),1),'om');
        set(w6,'XLim',[1,buffersize]);
        figure(2);
    end
    offset = offset + length(filtered_data(overlap1:end,:));
    pause
end

% add the last unprocessed portion

overlapbuffer = datasegment(size(datasegment,1)-overlap11: ...
    size(datasegment,1),:);
% buggy ?
% datasegment2 = [overlapbuffer;flipud(overlapbuffer)];

datasegment2 = [datasegment;flipud(overlapbuffer)];

[datasegment2,initvec] = adjovfl(datasegment2,initvec);

filtered_data = filter(firfiltb,1,datasegment2(:,channels));
filtered_data = abs(filtered_data);
[filtered_data2 ,filcond] = filter(avgfiltb,1,filtered_data,filcond);
filtered_data3 = sum(filtered_data2,2);

[dd,threshstate] = threshres(filtered_data3(overlap1+overlap/2-1,1), ...
    sthreshold,threshstate);
dd = (abs(dd) + offset) .* sign(dd);
sdetected = [sdetected;dd];

fclose(datafile);

% compensate the averaging filter delay
sdetected = (abs(sdetected) - avgfilterdelay) .* sign(sdetected);

%verdisp('done',verbose);

%%% sorter %%%%

%verdisp('sorting...',verbose);
minsp = floor(min_sw_period/1000 *  sampl);
maxsp = ceil(max_sw_period/1000 *  sampl);
minisp = ceil(min_isw_period/1000 *  sampl);

ss = sdetected(find(sdetected>0));
es = abs(sdetected(find(sdetected<0)));

ss = ss - regionexpand;
es = es + regionexpand;

swp = [];
lastswp = 0;
for s = ss'
    if lastswp < s
        e = find(es>=s+minsp & es<s+maxsp);
        if ~isempty(e)
            ee = e(1);
            swp = [swp;s,es(ee)];
            lastswp = es(ee);
        end
    end
end

swpn = [];

if (~isempty(swp))
    swptmp = swp(1,:);
    for ii=1:size(swp,1)-1
        swint = swp(ii+1,1)-swp(ii,2);
        if swint > minisp
            swpn = [swpn;swptmp];
            swptmp = swp(ii+1,:);
        else
            swptmp = [swptmp(1,1),swp(ii+1,2)];
        end
    end
    swpn = [swpn;swptmp];
    
    sm = floor(mean(swpn,2));
    sdetected = [swpn(:,1),sm,swpn(:,2)];
    %verdisp('done',verbose);
    if nargout<1
        swpsave(sdetected);
    end
else
    sdetected=[];
    
end
end

function [outbuffer,lastvec] = adjovfl(inbuffer,initvec)

% detect jumped buffer value due to overlap error (datamax)

[o1,o2] = find(abs(diff([initvec;inbuffer],1,1))>32000);

outbuffer = inbuffer;

% just incase ..
if (~isempty(o1))
    oo = sortrows([o1,o2],1);
    for ii = 1:length(o1)
        if oo(ii,1)==1
            outbuffer(1,oo(ii,2)) = initvec(oo(ii,2));
        else
            outbuffer(oo(ii,1),oo(ii,2)) = outbuffer(oo(ii,1)-1,oo(ii,2));
        end
    end
end

lastvec = inbuffer(end,:);
end

function [resvec,state] = threshres(data,thresh,state)

data = data(:);

tdata = [state;(data > thresh)];
dtdata = diff([tdata]);
resvec = find(dtdata).*sign(dtdata(find(dtdata)));
state = tdata(end);
end 