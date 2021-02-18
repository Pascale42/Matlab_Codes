% function DetectUPstates(FileBase, Name, channel, method, Save, Threshold)
%
% FileBase :: Root name of the experiment in a string
% Name :: Name of the structure where the UPstates will be detected
% channel :: Channel used for detection
% method :: 'detect' for regular detection
%                 'whiten' works well only with cortical LFP, 
% Save :: 0 or 1
% Threshold :: optional: by default 1.2 for 'whiten' and 0.1 for 'detect'
%
% Output::  UPstates structure
% UPstates.int :: beginning/end of UPstates in samples
% UPstates.len :: UPstates duration in msec
% UPstates.t :: middle of the UPstates in samples


function UPstates = DetectUPstates(FileBase, Name, channel, method, Save, Threshold)

if nargin < 6
    Threshold=1.2;
end

Par = LoadPar([FileBase '.xml']);
eeg=LoadBinary([FileBase '.eeg'], channel, Par.nChannels);

STA = load([FileBase '.sts.SWS']);
[eeg, orind]=SelectPeriods(eeg(:,:), STA, 'c', 1);


switch method
    
    case 'whiten'
        
        if nargin < 6
            Threshold=1.2;
        end
        
        eeg_new = WhitenSignal(eeg, Par.lfpSampleRate*2000,2);
        Std_new=std(eeg_new);
        
        Crossings = SchmittTrigger(eeg_new,Std_new*Threshold,Std_new*Threshold);
        
        CrossingsDiff=diff(Crossings);
        jump= CrossingsDiff ~= 1;
        CrossingsJumps=Crossings(jump);
        CrossingsJumpsDiff=diff(CrossingsJumps);
        
        Upfin=find(CrossingsJumpsDiff > 100); % end of UP-state
        Even=CrossingsJumps(Upfin);
        Even=[Even; CrossingsJumps(end)]; %  to get last point
        
        Updeb= Upfin+1; % beginning of UP-state
        Odd=CrossingsJumps(Updeb);
        Odd=[CrossingsJumps(1); Odd]; % to get first point
        
        
        % figure; plot(eeg); hold on; plot(eeg_new, 'm');scatter(Crossings, repmat(5e03, length(Crossings), 1), 'r')
        % scatter(IdUpdeb(UpGood), repmat(3e03, length(IdUpdeb(UpGood)), 1), 'g')
        % scatter(IdUpfin(UpGood), repmat(3e03, length(IdUpfin(UpGood)), 1), 'm')
        
        
    case 'detect'
        
        if nargin < 6
            Threshold=0.6;
        end
        
        % Get only slow oscillations
        feeg = ButFilter(eeg,2,[0.1 2]/(Par.lfpSampleRate/2),'bandpass');
        
        % Get windowed variance of the feeg
        Dt=100;
        Vf=zeros(length(feeg)- Dt,1);
        for t=1 : length(feeg) -Dt
            Vf(t) = var(feeg(t : t + Dt));
        end
        Vf=[zeros(Dt/2,1); Vf; zeros(Dt/2,1)];
        
        % find peaks in Vf
        [~,locs] = findpeaks(Vf, 'MINPEAKHEIGHT',Threshold*std(Vf));
        % figure; plot(Vf); hold on;
        % % offset values of peak heights for plotting
        % scatter(locs,pks+0.05,'k^','markerfacecolor',[1 0 0]);
        clear Vf
        
        % Get which transition
        Trans=ones(length(locs),1);
        for n=1:length(locs)
            
            if mean(feeg(locs(n)-Dt/2))-mean(feeg(locs(n)+Dt/2)) > 0
                Trans(n)=1;
                % U>D
            else
                Trans(n)=-1;
                % D>U
            end
        end
        
        % taking care of double detections:
        DoubleTrans=Trans(1:end-1)+Trans(2:end);
        badud=find(DoubleTrans == -2); % -2 = double D>U   and keep the last
        baddu=find(DoubleTrans == 2);  % 2 = double U>D
        baddu=baddu+1; % keep the first
        DoubleTrans= sort([baddu; badud]);
        id=~ismember([1:length(Trans)], DoubleTrans);
        Transitions=[locs(id) Trans(id)];
        clear id baddu badud DoubleTrans Trans locs n
        
        
        % start with beginning UPstate
        if Transitions(1,2) ~= 1
            Transitions=Transitions(2:end, :);
        end
        
        % end with end of UPstate
        if Transitions(end,2) ~= -1
            Transitions=Transitions(1:end-1, :);
        end
        
        Transitions=Transitions(:,1);
        
        Even = Transitions(2:2:length(Transitions));
        Odd = Transitions(1:2:length(Transitions));
        
end


%  minimal duration of  UP-state
len = Even - Odd; % duration all UP-state
Good = find(len >= 300); % get the ones long enough

% Ouput

% UPstates.int =[orind(Odd(Good)) orind(Even(Good))];     % !!!!!! ORIND !!!!!!
UPstates.int =[Odd(Good) Even(Good)];     % !!!!!! NO   ORIND !!!!!!

UPstates.len = len(Good) / Par.lfpSampleRate *1000; % in msec
UPstates.t = round((UPstates.int(:,1) + UPstates.int(:,2))/2); % middle of the burst
UPstates.par.thresh=Threshold;

UPstates.par.method=method;
UPstates.par.chan=channel;

if Save
    save([FileBase '.' mfilename '.' Name '.mat'], 'UPstates');
end




