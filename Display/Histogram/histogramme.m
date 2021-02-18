function out=histogramme(FileName,Threshold1,Threshold2,anatgrps,channels,anatgrpstest,...
    channelstest,channelstestfin,maxmean) %add anatgrpstest, channelstest
%maxmean=1->max
%maxmean=0->mean
    
B=[anatgrpstest,channelstest];
C=[channelstestfin];
if anatgrps==anatgrpstest & channels==channelstest
    error('Channel test = Channel');
end 

if isempty(C)       %si channelstestfin est vide
    if isempty(B)    %si anatgrpstest et channelstest vides
    FileBase=FileName;
    Par=LoadPar([FileBase '.xml']);
    Channels=Par.AnatGrps(anatgrps).Channels([channels])+1;
    %Channelstest=Par.AnatGrps(anatgrpstest).Channels([channelstest])+1;
    %EEG=LoadBinary([FileBase '.eeg'],Channels, Par.nChannels);
    out.Ripples=DetectRipples(FileBase, Channels,[],[100 200],[Threshold1 Threshold2]);

    
    else
    FileBase=FileName;
    Par=LoadPar([FileBase '.xml']);
    Channels=Par.AnatGrps(anatgrps).Channels([channels])+1;
    Channelstest=Par.AnatGrps(anatgrpstest).Channels([channelstest])+1;
    %EEG=LoadBinary([FileBase '.eeg'],Channels, Par.nChannels) ;
    out.Ripples=DetectRipples(FileBase,  Channels,Channelstest,[100 200],[Threshold1 Threshold2]);
   
    end
else
    for i=channelstest:1:channelstestfin

    FileBase=FileName;
    Par=LoadPar([FileBase '.xml']);
    Channels=Par.AnatGrps(anatgrps).Channels([channels])+1;
    Channelstest=Par.AnatGrps(anatgrpstest).Channels([i])+1;
    %EEG=LoadBinary([FileBase '.eeg'],Channels, Par.nChannels);
    out.Ripples=DetectRipples(FileBase,Channels,Channelstest,[100 200],[Threshold1 Threshold2]);
    
    out.nb(i-channelstest+1)=length(out.Ripples.t);
    
     if 1
      bar(out.nb);
     end 
    end
if maxmean 
        [m,j]=max(out.nb); %j=indice max
        Channelstest=Par.AnatGrps(anatgrpstest).Channels([j])+1;
        out.Ripples=DetectRipples(FileBase,Channels,Channelstest,[100 200],[Threshold1 Threshold2]);
        
    else
        moy=mean(out.nb);
        nb1=abs(out.nb-moy);
        [m, j]=min(nb1);       %indice plus proche de la moyenne
        Channelstest=Par.AnatGrps(anatgrpstest).Channels([j])+1;
        out.Ripples=DetectRipples(FileBase,Channels,Channelstest,[100 200],[Threshold1 Threshold2]);
       
        
    end
end
    out.FileBase=FileBase;
    out.Par=Par; 
    %out.Channels=Channels;
    %out.Channelstest=Channelstest;
    %out.EEG=EEG;
    out.Threshold=[Threshold1 Threshold2];
end 

% subplot(4,1,1);
% hist(ans1.zpow,one);en 
% subplot(4,1,2);
% hist(ans2.zpow,one);
% subplot(4,1,3);
% hist(ans3.zpow,one);
% subplot(4,1,4);
% hist(ans4.zpow,one);