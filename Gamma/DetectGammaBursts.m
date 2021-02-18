% 
% function [GammaBursts, GBduration, GBtime] = DetectGammaBursts(FileBase, State, Name, channel, FreqRange,Threshold)
%     
% FileBase :: Root name of the experiment in a string
% State :: 'THE, 'SWS'....
% Name :: Name of the structure where Gamma will be  detected
% channel :: Channel used for detection
% FreqRange :: to filter the LFP
% Threshold :: 3.5 by default
% 
% Output::     
% GammaBursts :: beginning/end of gamma bursts in samples
% GBduration :: gamma bursts duration in msec
% GBtime :: middle of the burst in samples


function [GammaBursts, GBduration, GBtime] = DetectGammaBursts(FileBase, State, Name, channel, FreqRange, Threshold)


Par=LoadPar([FileBase '.xml']);
% STA=load([FileBase '.sts.' State]); 
% eeg=LoadBinary([FileBase '.eeg'], channel, Par.nChannels);
% [eeg, orind]=SelectPeriods(eeg(:,:), STA, 'c', 1);
[eeg, orind]=GetEegState(FileBase, channel, 0, State);


feeg = ButFilter(eeg,2,FreqRange/(Par.lfpSampleRate/2),'bandpass'); 

weeg = WhitenSignal(eeg,Par.lfpSampleRate*2000,2);
fweeg = ButFilter(weeg,2,FreqRange/(Par.lfpSampleRate/2),'bandpass'); 
clear eeg;


if nargin < 6
Threshold = 2.5;
end

Stdev=std(feeg);
x=find(feeg  > Threshold*Stdev | feeg < -Threshold*Stdev); % thresholding 


%~~~~~~~~~~~DETECTION
x2=diff(x); 
jump= x2 ~= 1;
xj=x(jump);
xj2=diff(xj);

bfin=find(xj2 > 100); % end of gamma burst
xbf=xj(bfin); xbf=[xbf; xj(end)]; %  to get last point

bdeb= bfin+1; % beginning of gamma burst
xbd=xj(bdeb); xbd=[xj(1); xbd]; % to get first point

%  minimal duration of  gamma burst
bduree=xbf-xbd; % duration all burst
bgood=find(bduree > 100); % get the ones long enough

clear x x2 jump xj  bfin bdeb bduree

GammaBursts=[xbd(bgood) xbf(bgood)]; % beginning/end of gamma bursts
%~~~~~~~~~~~

Stdevw=std(fweeg);
x=find(fweeg  > Threshold*Stdevw | fweeg < -Threshold*Stdevw); % thresholding 



%~~~~WINDOWED VARIANCE

Dt=100;  %  window, best so far
V=zeros(length(eeg)- Dt,1);
for t=1 : length(eeg) -Dt
    V(t) = var(eeg(t : t + Dt));
end

% make V as long as fegg: add half window zeros each side
V=[zeros(Dt/2,1); V; zeros(Dt/2,1)];

% before making V the right length, had to adjust xbd and feeg:
% plot( feeg(Dt/2:length(feeg)-(Dt/2)) ); scatter(xbd(bgood)+Dt/2,
% repmat(2, length(xbd(bgood)),1))
%~~~~


figure
set(gca,'Position',[.02 .04 .96 .93])
plot(eeg,'k')
hold on
plot(V./max(V)*max(eeg) + 3000,'g')

plot(feeg +1500,'Color', [1.00	0.55	0.41])
scatter(xbd(bgood), repmat(-5*std(feeg), length(xbd(bgood)),1), 'og', 'LineWidth', 2)
scatter(xbf(bgood), repmat(-5*std(feeg), length(xbf(bgood)),1), 'or', 'LineWidth', 2) % threshold 2.5

scatter(xbd(bgood), repmat(-2*std(feeg), length(xbd(bgood)),1), '*g', 'LineWidth', 2)
scatter(xbf(bgood), repmat(-2*std(feeg), length(xbf(bgood)),1), '*r', 'LineWidth', 2)% threshold 2.

plot(fweeg .*3 -3000,'Color', [0.8,0.2,0.22])
scatter(xbd(bgood), repmat(-40*std(fweeg), length(xbd(bgood)),1), 'd', 'LineWidth', 2,'MarkerFaceColor',[0.55	0.53	0.31])
scatter(xbf(bgood), repmat(-40*std(fweeg), length(xbf(bgood)),1), 'd', 'LineWidth', 2,'MarkerFaceColor',[0.80	0.69	0.58])


% NOPE, do 2nd condition on the basis of variance
%
% Vnorm= V./max(V)*100 +1;
% feeg_corr= feeg(Dt/2:length(feeg)-(Dt/2) -1) ./Vnorm;
% figure; plot(feeg(Dt/2:length(feeg)-(Dt/2) -1));hold on; plot(feeg_corr,'k')
% 
% 
% Stdev=std(feeg_corr);
% x=find(feeg_corr  > Threshold*Stdev | feeg_corr < -Threshold*Stdev); % thresholding 





%~~~~~~~~~~~DETECTION
x2=diff(x); 
jump= x2 ~= 1;
xj=x(jump);
xj2=diff(xj);

bfin=find(xj2 > 100); % end of gamma burst
xbf=xj(bfin); xbf=[xbf; xj(end)]; %  to get last point

bdeb= bfin+1; % beginning of gamma burst
xbd=xj(bdeb); xbd=[xj(1); xbd]; % to get first point

%  minimal duration of  gamma burst
bduree=xbf-xbd; % duration all burst
bgood=find(bduree > 100); % get the ones long enough

clear x x2 jump xj xj2 bfin bdeb 
%~~~~~~~~~~~
GammaBursts=[xbd(bgood) xbf(bgood)]; % beginning/end of gamma bursts

plot(eeg,'k'); hold on; axis tight; plot(feeg +1500,'Color', [1.00	0.55	0.41])
scatter(xbd(bgood), repmat(5*std(feeg), length(xbd(bgood)),1), 'og', 'LineWidth', 2)
scatter(xbf(bgood), repmat(5*std(feeg), length(xbf(bgood)),1), 'or', 'LineWidth', 2)

plot(weeg)
plot(fweeg.*3 - 1000,'Color', [0.8,0.2,0.22])
scatter(xbd(bgood), repmat(-10*std(fweeg), length(xbd(bgood)),1), 'd', 'LineWidth', 2,'MarkerFaceColor',[0.55	0.53	0.31])
scatter(xbf(bgood), repmat(-10*std(fweeg), length(xbf(bgood)),1), 'd', 'LineWidth', 2,'MarkerFaceColor',[0.80	0.69	0.58])



GammaBursts=[xbd(bgood) xbf(bgood)]; % beginning/end of gamma bursts
GammaBursts=[orind(GammaBursts(:,1)) orind(GammaBursts(:,2))];
GBduration= (xbf(bgood) - xbd(bgood))/Par.lfpSampleRate *1000; % in msec
% WRONG !!!!! GBtime= round((GammaBursts(:,1) + GammaBursts(:,2))/2); % middle of the burst NO!!!
GBmiddle= round((GammaBursts(:,1) + GammaBursts(:,2))/2); % middle of the burst
GBtime= GammaBursts(:,1); % Beginning of the burst

MakeEvtFile(GBtime, [FileBase '.' Name '.' State '.gam.evt'],'gam',Par.lfpSampleRate,1);
save([FileBase '.' mfilename '.' Name '.' State '.mat'], 'GammaBursts', 'GBduration','GBtime','GBmiddle');

end
