%function [StTimes, StClu]=DetectArtifact(FileBase,ChannelToDet);
% or function [StTimes, StClu]=DetectArtifact(EEG, FileBase);
%detectes stimulus artifacts from .dat file
function [StTimes, StClu]= DetectArtifact(FileBase,ChannelToDet)



Plotwin  = [-10 25]; % msec
FreqRange = [1000 5000];
Par = LoadPar(FileBase);
Fs = Par.SampleRate;
histbuff = 100*Fs;
%EEG=readsinglech([FileBase '.dat'],Par.nChannels,ChannelToDet);
EEG=LoadBinary([FileBase '.dat'],ChannelToDet,Par.nChannels);

MyFilt=Sfir1(30, FreqRange/Fs*2);
FiltEEG=Filter0(MyFilt,EEG);
Sign = sign(abs(max(FiltEEG)) - abs(min(FiltEEG)));
figure(989)
clf
subplot(211)
hist(FiltEEG(1:histbuff),1000);
fprintf('click on the location of the amplitude threshold for peaks detection\n');
[DThr, y] = ginput(1);

if sign(DThr)~=Sign Sign=-Sign; end
PutativeTimes = LocalMinima(-Sign*FiltEEG,100,-Sign*DThr);

hist(FiltEEG(PutativeTimes),100);
isok =0;
while ~isok
    fprintf('click on the location of the final threshold \n');
    subplot(211)
    [FinalThr, y] = ginput(1);
    
    StTimes = PutativeTimes(find(-Sign*FiltEEG(PutativeTimes)<-Sign*FinalThr));
    fprintf('detected %d events\n',length(StTimes));
    PlotWin = [round(Plotwin(1)*Fs/1000):round(Plotwin(2)*Fs/1000)];
    seg = repmat(StTimes,1,length(PlotWin))+ repmat(PlotWin(:)',length(StTimes),1);
    datseg = EEG(seg);
    subplot(212)
    plot(PlotWin/Fs*1000,datseg,'Color','k');
    isok = ~input('do you want to redetect (0 or 1)?\n');
end        
ifcluster = input('Do you want to cluster the responses (1 -YES / 0 - NO) ?\n');
if ifcluster
    fprintf('click on the approximate region of the first EP (monosynaptic reponse)\n');
    fprintf('click first on left and the on right boundary\n');
    subplot(212)
    [EPboundary ,y ]= ginput(2); 
    EPboundary = EPboundary - Plotwin(1); 
    EPboundary = round(EPboundary*Fs/1000);
    if (EPboundary(2)<EPboundary(1)) EPboundary = fliplr(EPboundary);end
     EPwin = [EPboundary(1):EPboundary(2)];  
%    EPwin = EPwin + repmat(seg(:,1),1,length(EPwin));
    EPseg = datseg(:,EPwin);
    [EPAmp EPmin]= min(EPseg,[],2);
    

     subplot(211)
     hist(EPAmp,20);
     nclu = input('Tell me how many clusters you suspect ..\n');
     fprintf('click on their approximate locations');
     [clucenters,y] = ginput(3);
%     [pc, pv ] = pca(EPseg);
%     EPpr  = EPseg*pv(:,1:2);
%     plot(EPpr(:,1),EPpr(:,2),'.');
%     keyboard
    
    addpath('/u12/antsiro/matlab/toolboxes/netlab');
    mix = gmm(1,nclu,'spherical');   
    opt = zeros(15,1); opt(14)=100; opt(5)=1;
    mix.centres = clucenters;
    mix = gmminit(mix, EPAmp,opt);
    hold on
    Lines(mix.centres,[],'r');
    Lines(clucenters,[],'g');
    StClu = gmmclu(mix,EPAmp);
    AvEP = [];
    for i=1:nclu
        myst = find(StClu==i);
        AvEP(i,:) = mean(datseg(myst,:));
    end
    subplot(212)
    hold on
    plot(PlotWin/Fs*1000,AvEP,'LineWidth',2);
%    legend(num2str([1:nclu]'));
end    

if (nargout<1)
    if FileExists([FileBase '.stim'])
        system([ 'rm ' FileBase '.stim']);
        system([ 'rm ' FileBase '.stimclu']);
        system([ 'rm stim']);
    end
    vsave([FileBase '.stim'],StTimes);
    vsave([FileBase '.stimclu'],StClu);
    vsave('stim',StTimes);
    MakeEvtFile(StTimes,[FileBase '.stm.evt'],'stim',Fs,1);
end

% if nargout>0
%     StimTime = SpTimes;
% end