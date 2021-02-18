%function Out = PhaseOfTrain(Ph,Clu,Res,Periods,InOut, Brief,nClu)
% Out. Clui, pval, th0, R, phconfint, mu k, phhist, phbin
% trains empty ones will not be filled up
% NB! all output is in radians, degrees suck!
% Clu - cluster number from the Clu vector given as an input. This is handy if you run it on 
% the data that has gaps in clusters e.g. unique(Clu) = [1 3 4 6] 
function Out = PhaseOfTrain(Ph,Clu,varargin)
[Res,Periods,InOut,Brief,nClu] = DefaultArgs(varargin,{[],[],1,0,max(Clu)});

uClu = unique(Clu);
nClus = length(uClu);

if ~isempty(Res) && ~isempty(Periods)
	[Res, ind] = SelectPeriods(Res,Periods,'d',InOut);
	Ph=Ph(ind);Clu=Clu(ind);
end


Out = struct;
for myi=1:nClus
   % if ismember(myi,uClu)
        if sum(Clu==uClu(myi))<2
            continue;
        end
        myph = Ph(Clu==uClu(myi));
        Out(myi).Clu = uClu(myi);
        %myph = Ph(Clu==uClu(myi));
        [Out(myi).pval, th0, Out(myi).R] = RayleighTest(myph);
        Out(myi).th0 = mod(th0,2*pi);
        if ~Brief
            [cm, r] = CircConf(myph,0.01,200);
            Out(myi).phconfint = mod(r,2*pi)';
            [Out(myi).mu, Out(myi).k] = VonMisesFit(myph);
        end
       %myph = mod(myph,2*pi);
%        [Out(myi).phhist Out(myi).phbin] = hist(myph,[0:pi/18:2*pi]);
        [Out(myi).phhist, Out(myi).phbin] = hist([myph; myph+2*pi],[-pi:pi/18:3*pi]);
    %end
end

