%function TrigPower(filename, win,Resample)
function out =TrigPower(filename, varargin)

[win, Resample] = DefaultArgs(varargin, {300, 1});
winsample  = ceil(1250*win/1000/Resample);
nSamples=2*winsample+1;
load([filename '.osc.mat'],'pow');
pow = log(pow(1:Resample:end,:));
load([filename '.DetectPeaks.mat']);
peak = PeakPeak(filename);
sel = find(ismember(OutArgs(:,2),[1 2]));
peak(find(peak(:,2)==3),2)=2;
t = [OutArgs(sel,1);peak(:,1)];
g =[OutArgs(sel,2);peak(:,2)+2];
t= floor(t/Resample);
rem = load([filename '.rem']);
rem = round(rem/16/Resample);

[t, gi ] = SelectPeriods(t,rem,'d',0);
g=g(gi);

% so now we have 
labelp = {'delta', 'sp. spike', 'sp. peak', 'ripple peak'};

[avp stdp] = TriggeredAv(pow, winsample,winsample, t,g);

%[distrp t histp] = TriggeredDistr(pow, winsample,winsample, t,g,20,1, 20,0.005,1250);
time =linspace(-win, win,nSamples);
if (nargout>0)
    out{1} = avp;
    out{2} = stdp;
    out{3} = time;
    % out{4} =distrp;
    % out{5} = histp;
    
else
    figure(1)
%    pm = permute(cat(1, shiftdim(avp,-1), shiftdim(avp-stdp,-1), shiftdim(avp+stdp,-1)),[2 1 3 4]);
%    pm = SmoothMatrix(avp, 0.002);   
    PlotMatrix(time/1000, pm);
    ForAllSubplots('xlim([-19 19])');
    %     figure(2)
    %     ImageMatrix(time, histp, distrp);
end

