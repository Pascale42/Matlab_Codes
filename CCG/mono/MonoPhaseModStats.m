%function out = MonoPhaseMod(FileBase,fMode)

function out= MonoPhaseMod(FileBase,varargin)
[fMode] = DefaultArgs(varargin,{'compute'});
switch fMode
    case 'compute'
       % load([FileBase '.' mfilename '.mat']);
        
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
%        Randomiz = struct('nRand',1000,'Tau',5,'nSubs',50,'Alpha',[5 95]);
         Randomiz = struct('nRand',1000,'Tau',5,'nSubs',50,'Alpha',[5 95]);
        MaxLag = 10;
     
        load([FileBase '.thpar.mat'],'ThPh');
%         Ph = ThPh(round(Res/16));
%         [dummy ThInd] = histcI(Ph,linspace(-pi,pi,9));

        wcnt=1;
        for w=1:2
            if ~FileExists([FileBase '.sts.' States{w}]); continue; end
            Period =load([FileBase '.sts.' States{w}]);
            Period = Period*Par.SampleRate/Par.lfpSampleRate;

             [myRes, ind] = SelectPeriods(Res,Period,'d',1);
             myClu = Clu(ind);
             out(wcnt).monomod = MonoModSignif(ThPh,myRes,myClu,MaxLag,PairsToTest,Randomiz,Par.SampleRate);
            out(wcnt).Map = Map(:,2:3);
            out(wcnt).PairsToTest = PairsToTest;
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
        mout = load([FileBase '.MonoPhaseMod.mat']);
        nStates = length(out);
        figure(323232);
        for w=1:nStates
          nCells = length(out(w).monomod);
        
          for c=1:nCells
            clf
            subplot(221);cla
            tbin = mout.out(w).phccg.tbin(:,1);
            bar(tbin,mout.out(w).totccg.CCG(:,c));axis tight
            hold on
            plot(tbin,mout.out(w).totccg.smCCG(:,c),'LineWidth',3);axis tight
            plot(tbin,mout.out(w).totccg.AvShufCCG(:,c),'r'); 
            sigt = find(mout.out(w).totccg.PvalShufCCG(:,c)<0.01);
            ylim(ylim.*[1 1.1]);
            SigPer = tbin(Ind2Per(sigt,2)');
            %                SigPer = dotdot(SigPer,'+',[-1 1]*diff(tbin(1:2))/2);
            ShadeArea(SigPer,[]);
            xlim([-10 10]);
            title([num2str(mout.out(w).mono.From.ElClu(c,:)) ' -> ' num2str(mout.out(w).mono.To.ElClu(c,:))]);
            
            myout = out(w).monomod(c);
            if isempty(myout.lmodjit) | ~isstruct(myout.lmod) continue; end
            lag = myout.lmod.lag;
            pval = sum(dotdot(myout.lmodjit.r,'-',myout.lmod.r)>0,2)/1000;
            subplot(222)
            bar(lag,myout.lmod.r);axis tight
            subplot(224);
            bar(lag,-log10(pval));axis tight
            xlim([-10 10]);
            subplot(223)
            nx = size(myout.lmod.mhst,1);
            mhst = NaN*zeros(nx+3,8);
            mhst(1:nx/2,:) = myout.lmod.mhst(1:nx/2,:);
            mhst(end-nx/2+1:end,:) = myout.lmod.mhst(nx/2+1:end,:);
            mhst(:,[1 8]) = repmat(sum(mhst(:,[1 8]),2),1,2);
            nlag = lag(1):lag(end);
            imagesc(lag,linspace(-pi,pi,8),unity(mhst'));
            [x y b] = ginput(1);
            if b==2
                return;
            end
          end
 
        end
        
        

end
