%function SpectralCellAnal(FileBase,ElClu, States, FreqRange,s2f)
% plots the spectral analysis for a cell (uses s2f file)
% you have to run Spikes2Field beforehand (may take quite some time)
% you can narrow plot to particular states
function SpectralCellAnal(varargin)
[FileBase, ElClu, States,FreqRange,s2f] = DefaultArgs(varargin,{[],[1 2],[],[1 20],[]});
if isempty(FileBase)
    [dummy FileBase] = fileparts(pwd);
end
if isempty(s2f)
    load([FileBase '.s2f'],'-MAT');
end
Par =LoadPar([FileBase '.par']);
EegCh = load([FileBase '.eegseg.par']);
nEegCh = length(EegCh);

tot = CatStruct(s2f,{'State','ElClu'});
if isempty(States)
    States = tot.State;
end
%^States = intersect(States,tot.State);

if isstr(ElClu)
    El = find(strcmp(Par.ElecLoc,ElClu));
    myi = find(ismember(sq(tot.ElClu(:,1,1)),El));
    ElClu = tot.ElClu(myi,:,1);
end
nCells = size(ElClu,1);
hlns = [];

nStates = length(s2f);
nFreq = size(FreqRange,1);
figure(666); clf

for j=1:nCells

    for i=1:nStates
        %^i = find(strcmp({'REM','SWS','RUN'},States{i}));
        if isempty(s2f(i).State)
            continue;
        end
        switch s2f(i).State
            case 'REM'
                mycol = 'g';
            case 'SWS'
                mycol = 'b';
            case 'RUN'
                mycol = 'r';
        end
        
        %in which area is the cell?
        CellLoc = Par.ElecLoc{ElClu(j,1)};
        CellInd = find(s2f(i).ElClu(:,1)==ElClu(j,1) & s2f(i).ElClu(:,2)==ElClu(j,2));
        switch CellLoc
            case 'h'
                RefCh = 1; 
            case 'c'
                RefCh = 2;
        end
        %plot power
        for f=1:nFreq
            subplot(3,nFreq, f)

            fi = find(s2f(i).f>FreqRange(f,1) & s2f(i).f<FreqRange(f,2));

            plot(s2f(i).f(fi), sq(s2f(i).y(fi,nEegCh+CellInd,nEegCh+CellInd)),mycol);
            hold on
            title([FileBase ' ' num2str(ElClu(j,:))],'FontSize',4);
            axis tight
            grid on
            ylabel('Unit power');

            %plot coherence to hpc and cx
            subplot(3,nFreq, nFreq+f)
            hlns(end+1) = plot(s2f(i).f(fi), sq(s2f(i).y(fi,RefCh,nEegCh+CellInd)),mycol);
            hold on
            %        plot(s2f(i).f, sq(s2f(i).y(:,2,nEegCh+CellInd)),[mycol '-.']);

            plot(s2f(i).f(fi), abs(sq(s2f(i).yerr(fi,RefCh,nEegCh+CellInd,1))),[mycol '--']);
            grid on
            %title([FileBase ' El=' num2str(ElClu(j,1))  ,' ,Clu='  num2str(ElClu(j,2)) ' (' CellLoc ')']);
            ylabel('Coherence ');
            axis tight

            %now phase
            subplot(3,nFreq, 2*nFreq+f)

            %find where coherence is above the
            bf = find(sq(s2f(i).y(:,RefCh,nEegCh+CellInd)) <= abs(sq(s2f(i).yerr(:,RefCh,nEegCh+CellInd,1))));
            phi = mod(unwrap(sq(s2f(i).phi(:,RefCh,nEegCh+CellInd))),2*pi);
            phierr= sq(s2f(i).phierr(:,RefCh,nEegCh+CellInd));
            phi(bf)=NaN; phierr(bf)=NaN;
            errorbar(s2f(i).f(fi), unwrap(phi(fi)), phierr(fi),mycol); hold on
            %        errorbar(s2f(i).f(gf), mod(sq(s2f(i).phi(gf,RefCh,nEegCh+CellInd)),2*pi), sq(s2f(i).phierr(gf,RefCh,nEegCh+CellInd)),[mycol '.']); hold on
            axis tight
            grid on
            ylabel('phase (Rad)');
            xlabel('Frequencey (Hz)');

        end %freq loop
    end %states loop
    %    ForAllSubplots(sprintf('xlim([ %d %d])',FreqRange(1),FreqRange(2)));
    [x y b] = PointInput(1);
    if b==3
        reportfig(gcf,'SpectralCellExamples',0,sprintf('File %s, El=%d, Clu=%d',FileBase, ElClu(j,1), ElClu(j,2)));
    elseif b==2
        return;
    end
    clf
end

% subplot(3,nFreq,1)
% legend(hlns,tot.State);
