%function out = MonoPhaseMod(FileBase,fMode)

function out= MonoPhaseMod(FileBase,varargin)
[fMode] = DefaultArgs(varargin,{'compute'});
switch fMode
    case 'compute'
        load([FileBase '.' mfilename '.mat']);
        
        if ~FileExists([FileBase '.mono-c'])
            fprintf('no mono pairs\n');
            out = [];
            return
        end
        States = {'REM','RUN'};
        load([FileBase '.mono-c'],'-MAT');
        if isempty(mono.From.ElClu)
            out = [];
            fprintf('no mono pairs\n');
            return
        end

       
        Par = LoadPar([FileBase '.xml']);
        El = find(strcmp(Par.ElecLoc,'c'));

        % phmod = load([FileBase '.PhaseFrAmpAnal.mat']);
        % phmodElClu = CatStruct(phmod.OutArgs,{'El','Clu'});

        %load([FileBase '.s2s'],'-MAT');
        %load([FileBase '.NeuronQuality.mat']);
        %OutArgs = CatStruct(OutArgs);
        %mycnq = find(ismember(OutArgs.ElNum,El));
        %nq = SubsetStruct(OutArgs, mycnq,1);

        %[el clu spkwdths ie] = textread([FileBase '.type-c'],'%d %d %s %d');

        n=size(mono.From.ElClu,1);
        mutind = [];
        muted = zeros(n,1);
        for ii=1:n
            if mono.To.Type==0 | muted(ii)==1
                continue;
            else
                tmp = find(mono.From.ElClu(:,1)==mono.To.ElClu(ii,1) & mono.From.ElClu(:,2)==mono.To.ElClu(ii,2) ...
                    & mono.To.ElClu(:,1)==mono.From.ElClu(ii,1) & mono.To.ElClu(:,2)==mono.From.ElClu(ii,2));
                if ~isempty(tmp)
                    mutind(end+1) = tmp;
                    muted(tmp) = 1;
                end
            end
        end
        mono.From.ElClu(mutind,:)=[];
        mono.To.ElClu(mutind,:)=[];
        mono.From.Type(mutind)=[];
        mono.To.Type(mutind)=[];

        n=size(mono.From.ElClu,1);

        [myElClu ui uj]= unique([mono.From.ElClu ; mono.To.ElClu],'rows');
        myEl = unique(myElClu(:,1));
        [Res,Clu,Map ] =  LoadCluRes(FileBase,myEl,myElClu);

        for ii=1:n
            Pre = find(Map(:,2) == mono.From.ElClu(ii,1) & Map(:,3) == mono.From.ElClu(ii,2));
            Post = find(Map(:,2) == mono.To.ElClu(ii,1) & Map(:,3) == mono.To.ElClu(ii,2));
            PairsToTest(ii,:)= [Pre Post];
        end

        %keyboard
        Jitter = 5*Par.SampleRate/1000;
        Randomiz = struct('Type','jitter','nRand',1000,'Tau',Jitter,'Alpha',[5 95]);
        BinSize = 10;
        HalfBins = round(20*Par.SampleRate/1000/BinSize); % to have 20msec total

        load([FileBase '.thpar.mat'],'ThPh');
        Ph = ThPh(round(Res/16));
        [dummy ThInd] = histcI(Ph,linspace(-pi,pi,9));

        % keyboard
%         for ii=1:n
%             myRes1 = SelectPeriods(Res(Clu==PairsToTest(ii,1)),Period,'d',1);
%             myRes2 = SelectPeriods(Res(Clu==PairsToTest(ii,2)),Period,'d',1);
%             %myPh = ThPh(round(myRes1/16));
%             [T,G] = CatTrains({myRes1,myRes2},
%             maxlag = 3*Par.SampleRate;
%             [xx xi yi] = NearestNeighbour(x,y,side, mindist);
%         
%         end%loop htough pairs
        wcnt=1;
        for w=1:2
            if ~FileExists([FileBase '.sts.' States{w}]); continue; end
            Period =load([FileBase '.sts.' States{w}]);
            Period = Period*Par.SampleRate/Par.lfpSampleRate;

            out(wcnt).totccg = CCGSignif(Res,Clu,BinSize,HalfBins,Par.SampleRate,'count',PairsToTest,Randomiz,10,Period);

           % out(wcnt).monomod =  MonoModSignif(ThPh,Res,Clu,BinSize,HalfBins,Par.SampleRate,'count',PairsToTest,Randomiz,5,Period);
            
%             for jj=1:max(ThInd)
%                 mRes= Res(ThInd==jj);
%                 mClu= Clu(ThInd==jj);
%                 out(wcnt).phccg(jj) = CCGSignif(mRes,mClu,BinSize,HalfBins,Par.SampleRate,'count',PairsToTest,Randomiz,5,Period);
% 
%             end
            [dummy, ind] = SelectPeriods(Res,Period,'d',1);
            out(wcnt).phhist = accumarray([ThInd(ind),Clu(ind)],1);
            out(wcnt).phccg = CatStruct(out(wcnt).phccg);
            out(wcnt).Map = Map(:,2:3);
            out(wcnt).PairsToTest = PairsToTest;
            out(wcnt).Ranodmiz = Randomiz;
            out(wcnt).State = States{w};
            out(wcnt).mono = mono;
            wcnt=wcnt+1;
        end

        save([FileBase '.' mfilename '.mat'], 'out');
    case 'display'
        if ~FileExists([FileBase '.' mfilename '.mat'])
            fprintf('no mono pairs here\n');
            return
        end
        load([FileBase '.' mfilename '.mat']);
        nStates = length(out);
        for w=1:nStates
          %  fprintf('State %s\n',out(w).State);
            tbin = out(w).phccg.tbin(:,1);
            if ndims(out(w).phccg.ZScore)==3
                nPairs = size(out(w).phccg.CCG,2);
                nphbin = size(out(w).phccg.CCG,3);
            else
                nPairs =1;
                nphbin = size(out(w).phccg.CCG,2);
            end
            phbin = linspace(-pi,3*pi,nphbin*2)*180/pi;
            myt = find(tbin>-10 & tbin<10);
            thhist = out(w).phhist;
            thhist = dotdot(thhist,'/',sum(thhist));
            ny= 2; nx =3;
            figure(3232)
            for ii=1:nPairs
                clf
                %which cells are those:
                myclu(1) = find(out(w).Map(:,1)==out(w).mono.From.ElClu(ii,1) & out(w).Map(:,2)==out(w).mono.From.ElClu(ii,2));
                myclu(2) = find(out(w).Map(:,1)==out(w).mono.To.ElClu(ii,1) & out(w).Map(:,2)==out(w).mono.To.ElClu(ii,2));
%                subplotfit(ii,nPairs);
                %%%%%% plot the CCG between pair members
                subplot(ny,nx,1);
                bar(tbin,out(w).totccg.CCG(:,ii));axis tight
                hold on
                plot(tbin,out(w).totccg.smCCG(:,ii),'LineWidth',3);axis tight
                plot(tbin,out(w).totccg.AvShufCCG(:,ii),'r');
%                plot(tbim,out(w).totccg.AvShufCCG(:,ii)+out(w).totccg.StdShuffCCG(:,ii)*3,'r--');
 %               plot(tbim,out(w).totccg.AvShufCCG(:,ii)-out(w).totccg.StdShuffCCG(:,ii)*3,'r--');
                sigt = find(out(w).totccg.PvalShufCCG(:,ii)<0.01);
                ylim(ylim.*[1 1.1]);
                SigPer = tbin(Ind2Per(sigt,2)');
%                SigPer = dotdot(SigPer,'+',[-1 1]*diff(tbin(1:2))/2);
                ShadeArea(SigPer,[]);
                xlim([-10 10]);
                title([num2str(out(w).mono.From.ElClu(ii,:)) ' -> ' num2str(out(w).mono.To.ElClu(ii,:))]);

                flds = {'','CCG','smCCG','AvShufCCG','ZScore'};
                for k=2:length(flds)
                    subplot(ny,nx,k);
                    if nPairs>1
                        mat = sq(out(w).phccg.(flds{k})(myt,ii,:));
                    else
                        mat = sq(out(w).phccg.(flds{k})(myt,:));
                    end
                    if k<5
                        mat = dotdot(mat,'/',sum(mat));
                    end
                    imagesc(tbin(myt),phbin,[mat mat]');
%                    imagesc(tbin(myt),phbin,conv2(gausswin(3),gausswin(3),[mat mat]','same'));
                    axis xy
                   set(gca,'YTick',[-180:180:540]);
                    title(flds{k});
                end
                subplot(ny,nx,6);
                stairs(repmat(thhist(:,myclu(1)),[2 1]),phbin,'b');
                hold on
                stairs(repmat(thhist(:,myclu(2)),[2 1]),phbin,'r');
                %axis tight
                ylim([-180 540]);
                set(gca,'YTick',[-180:180:540]);
           %     keyboard
                %waitforbuttonpress
                [x y b] = ginput(1);
                if b==3
                    keyboard
                end
            end
%            ForAllSubplots('xlim([-10 10])')
%            ForAllSubplots('caxis([-3 5])')

            
        end
        
        

end
