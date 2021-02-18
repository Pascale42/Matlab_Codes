%function [TrLag TrInd TrClu] = TrigRasters(Trig, Window, Res,Clu, SampleRate, CluSubplots, UseTrigOrder)
% plots the scatter plot of spikes in Res (colored by Clu)
% relative to the trigger Trig within a Window (in samples)
% two flags control the behavior of the plot:
% CluSubplots - if 1 it will make a separate subplot for each cluster,
% otherwise will put all in one coloring the dots of each cluster
% UseTrigOrder - if 0 Y axis will be the time of trigger occurrence,
% if 1 - the order index that triggers have. You can also provide it as the
% vector of length(Trig) and then sorting will be done according to the
% vector values
% this is usefull if you want to see how some aspect of the trigger (e.g.
% if its a ripple, then ripple power) has on the units behavior, then you
% presort the triggers by your favorite variable before calling this function
% REMINDER: make a cumulative CCG based on the CumulatingIndex vector
% provided for the reference (e.g. strength of event, gamma power, position
% etc ..if it is time - becomes cumulative JPSTH)
function [TrLag TrInd TrClu] = TrigRasters(Trig, Window, Res,varargin)

[Clu,SampleRate,CluSubplots,UseTrigOrder] = DefaultArgs(varargin,{ones(length(Res),1),1250,1,0});

nTrig = length(Trig);
%nClu = max(Clu);

if length(UseTrigOrder)>1
    %then it is the variable according to which we need to sort triggers
    if length(Trig)==length(UseTrigOrder)
        [dummy TrigSortIndOrig] = sort(UseTrigOrder);
        [dummy TrigSortInd] = sort(TrigSortIndOrig);
        UseTrigOrder=1;
    else
        error('length of Trig not equal to UseTrigOrder');
    end
else
    [Trig TrigSortInd] = sort(Trig);
end

%[T,G] = CatTrains({Trig, Res},{1,Clu});

T = [Trig(:); Res(:)]; 
G = [ones(length(Trig),1); Clu+1];

%[ccg tbin pairs] = CCG(T,G,Window,0,SampleRate,[1:max(G)]);
nClu = max(Clu);
Pairs2Use = [ones(nClu,1) [2:nClu+1]'];
[ccg tbin pairs] = CCG_part(T,G,Window,0,SampleRate,Pairs2Use);

TrPairs = pairs(G(pairs(:,1))==1 & G(pairs(:,2))>1,:);

%time lag from the trigger to the spike
TrLag = diff(T(TrPairs),1,2)/SampleRate*1000;

TrClu = G(TrPairs(:,2))-1;
TrTime = T(TrPairs(:,1));

%NB: Trig has to be sorted before using ismember and unique!
if UseTrigOrder
    [uTrTime dummy uTrInd] = unique(TrTime);
    TrInd = find(ismember(Trig,uTrTime));
    TrInd = TrigSortInd(TrInd(uTrInd));
else
     if nargout<2
         TrInd = TrTime./SampleRate*1000;
     end
%     TrInd =[1:length(TrTime)]';
      TrInd = TrPairs(:,1);
      
end


if nargout<2

    [uTrClu dummy TrClu] = unique(TrClu);
    nClu = max(uTrClu);
    if CluSubplots<2
        figure;clf
    end
    cols = colormap;
    cols = cols([1:round(64/nClu):64],:);
    for i=1:nClu
        myClu = find(TrClu==i);
        if ~isempty(myClu)
            if CluSubplots==1
                subplotfit(i,nClu);
                col = 'k';
            elseif CluSubplots==0
                hold on
                col = cols(i,:);
            elseif CluSubplots==2
                col = 'k';
            end
            h=plot(TrLag(myClu), TrInd(myClu),'.','Color',col,'MarkerSize',5);
            axis tight
            if CluSubplots==1
                title(num2str(['Cell # ' num2str(uTrClu(i))]));
            end
        end
    end
    if ~CluSubplots
        legend(num2str(uTrClu));
    end
    TrLag = h;
else
%    TrInd = TrigSortIndOrig(TrInd);
end


