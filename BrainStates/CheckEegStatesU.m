function CheckEegStates(FileBase,varargin)
%function CheckEegStates(FileBase,State)
State = DefaultArgs(varargin,{[]});

% auxil. struct for gui
global gCheckEegStates
gCheckEegStates = struct;

%constants
UseRulesBefore = 0; % flag to switch heuristics for periods cleaning after automatic segmentation
MinLen=5; %seconds
eSampleRate = 1250; % eeg sampling rate ..assume default
Par = LoadPar([FileBase '.par']);
nChannels = Par.nChannels;
Window = 1; %sec - for spectrogram calculation
FreqRange = [1 100]; % Hz - range of freq. for the spectrogram

if ~isempty(State)
    % load segmentation results
    if FileExists([FileBase '.sts.' State])
        Per = load([FileBase '.sts.' State]);
    else
        Per = SelectStates(FileBase,State,eSampleRate*MinLen);
    end       
    if isempty(Per)
        fprintf('No segments in %s state',State)
        %return;
        Per = [];
    end
else
    Per = [];
    State = '';
end
if FileExists([FileBase '.eegsegU.mat'])
    load([FileBase '.eegsegU.mat']); % load whitened spectrograms from EegSegmentation results
    Channels = load([FileBase '.eegseg.par'])+1; % loading channels to display as traces
else
    ch = inputdlg({'Enter channels to use'},'Channels',1,{'1'});
    Channels = num2str(ch{1});
    % now compute the spectrogram
    Eeg = readmulti([FileBase '.eeg'], nChannels, Channels);
    nFFT = 2^round(log2(2^11)); %compute nFFT according to different sampling rates
    SpecWindow = 2^round(log2(Window*eSampleRate));% choose window length as power of two
    weeg = WhitenSignal(Eeg,[],1);
    [y,f,t]=mtcsglong(weeg,nFFT,eSampleRate,SpecWindow,[],2,'linear',[],FreqRange);

end

t = (t(2)-t(1))/2 +t;

% computer the/del ratio and detect transitions automatically - not used at
% the momnet, maybe later
[thratio] = TDRatioAuto(y,f,t,MinLen);
%[thratio, ThePeriods] = TDRatioAuto(y,f,t,MinLen);

%now apply the rules to filter out junk states or make continuous periods
% to be implemented later
if UseRulesBefore
    switch State
        case 'REM'

    end
end

% fill the global structure for future use
gCheckEegStates.Channels = Channels;
gCheckEegStates.nChannels = nChannels;
gCheckEegStates.FileBase = FileBase;
gCheckEegStates.State = State;
gCheckEegStates.t = 100; %is seconds
gCheckEegStates.trange = [t(1) t(end)];
gCheckEegStates.Periods = Per/1250; % in seconds
gCheckEegStates.Mode = 't';
gCheckEegStates.nPlots=length(Channels)+2;%number of plots is fixed to 3: 2 spectrograms and one zoomin trace. If changed reuiqres fixes in _aux
gCheckEegStates.lh=cell(gCheckEegStates.nPlots,1);
gCheckEegStates.Window = Window*eSampleRate*2;
gCheckEegStates.SelLine=[];
gCheckEegStates.cposh=cell(gCheckEegStates.nPlots,1);
gCheckEegStates.FreqRange = FreqRange;
gCheckEegStates.newl=[];
gCheckEegStates.tstep = t(2)-t(1);
gCheckEegStates.coolln = [];
gCheckEegStates.LastBut = 'normal';

% create and configure the figure
gCheckEegStates.figh = figure('ToolBar','none');
set(gCheckEegStates.figh, 'Position', [3 828 1276 620]); %change Postion of figure if you have low resolution
set(gCheckEegStates.figh, 'Name', 'CheckEegStates : traces');
set(gCheckEegStates.figh, 'NumberTitle', 'off');


% put the uitoolbar and uimenu definitions here .. may require rewriting
% some callbacks as actions rather then cases of actions (e.g. key
% pressing)


%now do the plots

for ii=1:gCheckEegStates.nPlots-2
    subplot(gCheckEegStates.nPlots,1,ii);
    imagesc(t,f,log(sq(y(:,:,ii)))');axis xy; ylim([1 10]);
    hold on
    if ii==1
        title('Spectrogram'); ylabel('Frequency (Hz)');
    end
end
%  subplot(gCheckEegStates.nPlots,1,2)
%  plot(t,thratio);axis tight;
%  set(gca,'YTick',[]);
%  hold on
%  ylabel('Theta/Delta raio'); xlabel('Seconds');

subplot(gCheckEegStates.nPlots,1,gCheckEegStates.nPlots-1)
myt = [1:length(rippow)]/125;
plot(myt,rippow(:,1));axis tight

subplot(gCheckEegStates.nPlots,1,gCheckEegStates.nPlots)
CheckEegStatesU_aux('traces'); % plot the current eeg traces
hold on

%plots lines
if ~isempty(Per)
CheckEegStatesU_aux('lines');
end
% assign functions for mouse and keyboard click
set(gCheckEegStates.figh,'WindowButtonDownFcn','CheckEegStatesU_aux(''mouseclick'')');
set(gCheckEegStates.figh,'KeyPressFcn', 'CheckEegStatesU_aux(''keyboard'')');

%msave([FileBase '.states.' State],round(ThePeriods*eSampleRate));
return




function [thratio, ThePeriods] = TDRatioAuto(y,f,t,MinLen)

%automatic theta periods detection just using the thetaratio
thfin = find(f>6 & f<9);
thfout = find(f<5 | (f> 12& f<25));
thratio = log(mean(sq(y(:,thfin,1)),2))-log(mean(sq(y(:,thfout,1)),2));

if nargout>1
    nStates =2;
    % fit gaussian mixture and to HMM - experimental version .. uses only thetaratio
    [TheState thhmm thdec] = gausshmm(thratio,nStates,1,0);

    for i=1:nStates
        thratio_st(i) = mean(thratio(TheState==i));
    end

    [dummy TheInd] = max(thratio_st);
    InTh = (TheState==TheInd);
    DeltaT = t(2)-t(1);
    MinSeg = round(MinLen/DeltaT);

    TransTime = ThreshCross(InTh,0.5,MinSeg);
    ThePeriods = t(TransTime);
end
