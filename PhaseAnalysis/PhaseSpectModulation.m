function out = PhaseSpectModulation(FileBase,varargin)
%function PhaseSpectModulation(FileBase,fMode, ElLoc, FreqRange,State,FrThreshold,Window)
[fMode, ElLoc, FreqRange,States,FrThreshold,Window] = DefaultArgs(varargin,{'compute','c',[30 150], {'REM','RUN'},5,2^7});

Par = LoadPar([FileBase '.xml']);
MinPeriod = 5; %seconds for theta periods selection
if strcmp(ElLoc,'all')
    El = [1:Par.nElecGps]; ElClu = [];
elseif isstr(ElLoc(1))
    El = find(strcmp(Par.ElecLoc,ElLoc));ElClu = [];
else
    El = unique(ElLoc(:,1));ElClu = ElLoc;
end
if ischar(States); States = {States}; end;
nSt = length(States);

[RepCh ifTetr Info] = RepresentChan(FileBase);
nRepCh = length(RepCh);
eFs = 1250;

switch fMode
    case 'compute'
        [AllRes,AllClu,Map] = LoadCluRes(FileBase,El,ElClu);
        if isempty(AllRes)
            fprintf('no units in %s in this file \n',ElLoc);
            return;
        end
        AllRes = round(AllRes/16)+1;
     
        load([FileBase '.thpar.mat']);
        %loop through States
        for s=1:nSt
            Period = SelectStates(FileBase,States{s},MinPeriod);
            if isempty(Period)
                fprintf('No state %s\n',States{s});
                continue;
            end
            PerLen = round(sum(diff(Period,1,2)/1250));
            fprintf('Processing state %s of duration %d seconds\n',States{s},PerLen);

            [Res ind]=SelPerDiscr(AllRes,Period,1,1);
            Clu=AllClu(ind);
            [uclu dummy Clu] = unique(Clu);
            myMap = Map(uclu,:);
            uClu = unique(Clu);
            eeg = LoadBinary([FileBase '.eeg'],RepCh,Par.nChannels,Period)';
            myThPh = SelectPeriods(ThPh,Period,'c',1);
            fr = FiringRate(Res,Clu,[],1250);
            gfr = find(fr>FrThreshold);
            phbin = linspace(-pi,pi,6);
            nbins = length(phbin)-1;
            [hh PhInd] = histcI(myThPh(Res),phbin);
            coh=[]; ph=[];powf=[];powu=[];
            %            h = waitbar(0,'Please wait...');
            %                waitbar(ii/nRepCh,h);
            for jj=1:length(gfr);
                myRes = Res(Clu==gfr(jj));
                myPhClu = PhInd(Clu==gfr(jj));
                badph = myPhClu==0;
                myRes(badph)=[];
                myPhClu(badph)=[];
                myuClu = unique(myPhClu);
                if length(myRes)==0 continue; end
                if 1
                    for ii=1:nRepCh
                        [y, out(s).f,  phi, yerr, phierr, phloc] = mtptchd(eeg(:,ii),myRes,myPhClu,...
                            2^9,1250,2^nextpow2(Window),[],6,[],[],FreqRange,[],[],[],2);
                        out(s).powf(:,ii) = sq(y(:,1,1));
                        for k=1:length(myuClu)
                            out(s).powu(:,myuClu(k),jj) = sq(y(:,1+k,1+k));
                        end
                        out(s).coh(:,myuClu,ii,jj) = sq(y(:,1,2:end));
                        out(s).ph(:,myuClu,ii,jj) = sq(phi(:,1,2:end));
                        out(s).plph(:,myuClu,ii,jj) = sq(angle(phloc(:,1,2:end)));
                        out(s).ploc(:,myuClu,ii,jj) = sq(abs(phloc(:,1,2:end)));
                        out(s).yerr(:,myuClu,ii,jj) = sq(yerr(:,1,2:end,1));
                    end %channels
                
                else
                    %do all channels at once
                        [y, out(s).f,  phi, yerr, phierr, phloc] = mtptchd(eeg,myRes,myPhClu,...
                            2^9,1250,2^nextpow2(Window),[],6,[],[],FreqRange,[],[],[],2);
                        for k=1:nRepCh
                            out(s).powf(:,k) = sq(y(:,k,k));
                        end
                        for k=1:length(myuClu)
                            out(s).powu(:,myuClu(k),jj) = sq(y(:,nRepCh+k,nRepCh+k));
                        end
                        out(s).coh(:,myuClu,:,jj) = sq(y(:,1:nRepCh,nRepCh+1:end));
                        out(s).ph(:,myuClu,:,jj) = sq(phi(:,1:nRepCh,nRepCh+1:end));
                        out(s).plph(:,myuClu,:,jj) = sq(angle(phloc(:,1:nRepCh,nRepCh+1:end)));
                        out(s).ploc(:,myuClu,:,jj) = sq(abs(phloc(:,1:nRepCh,nRepCh+1:end)));
                        out(s).yerr(:,myuClu,:,jj) = sq(yerr(:,1:nRepCh,nRepCh+1:end,1));
                    
                end
            end %cells
            %           close(h);
            out(s).ElClu = myMap(gfr,2:3);
            out(s).FirRate = fr(gfr);
            out(s).State = States{s};
        end %loop through states

        if nargout<1
            OutArgs = out;
            save([FileBase '.' mfilename '.' States{s} '.mat'],'OutArgs');
        end


    case 'display'
        load([FileBase '.' mfilename '.mat']);
        if isempty(OutArgs(1).f) OutArgs(1)=[]; end
        phbin = linspace(-pi,pi,6);
        f = OutArgs(1).f; myf = find(f>30 & f<100);
        ni=2; nj=2;
        chi = 1;
        for s=1:length(OutArgs)
            nCells = size(OutArgs(s).ElClu,1);
            for i=1:nCells

                if OutArgs(s).FirRate(i)<5 continue; end
                while 1
                    figure(2126);
                    clf
                    subplot(ni,nj, 1)
                    %        npow = log(sq(OutArgs.y(:,myf,1,1)));
                    %npow = sq(OutArgs.y(:,myf,1,1)).*repmat(f(myf).^2',length(phbin),1);
                    npow = log10(sq(OutArgs(s).powf));
                    imagesc(f(myf),[1:nRepCh],npow');axis xy
                    title('power of eeg');

                    subplot(ni,nj, 2)
                    %npowspk = sq(OutArgs(s).y(:,myf,2,2)) ./repmat(interp1(OutArgs(s).tbin,sq(OutArgs(s).ccg(:,1,2)), OutArgs(s).phbin),length(myf),1)';
                    npowspk = sq(OutArgs(s).powu(:,:,i));
                    imagesc(phbin,f(myf),npowspk);axis xy
                    title('power of spikes');

                    % coh: nf x nph x nch x ncells
                    subplot(ni,nj, 3)
x``                    imagesc(phbin,f(myf),(sq(OutArgs(s).coh(myf,:,chi,i))));axis xy
                    title('coherence spks to eeg');
                    hold on
                    %phvar = 1./(sq(OutArgs(s).phierr(:,myf(1:5:end),1,2))+eps);
%                     u = real(exp(sqrt(-1)*sq(OutArgs(s).ph(myf(1:5:end),:,chi,i))));
%                     v = imag(exp(sqrt(-1)*sq(OutArgs(s).ph(myf(1:5:end),:,chi,i))));
%                     quiver(phbin,f(myf(1:5:end)),u,v,1, 'Color','k');

                    %         subplot(ni,nj, 4)
                    %         bar(OutArgs(s).tbin, OutArgs(s).ccg(:,1,2));axis tight
                    %         yl = ylim;    ylim([mean(yl) 1.3*yl(2)]);
                    %         title('spike ccg');

                    %plot again with normalization
                    subplot(ni,nj, 4)
                    imagesc(phbin,f(myf),(sq(OutArgs(s).ploc(myf,:,chi,i))));axis xy
                    title('coherence spks to eeg');
                    hold on
                    %phvar = 1./(sq(OutArgs(s).phierr(:,myf(1:5:end),1,2))+eps);
%                     u = real(exp(sqrt(-1)*sq(OutArgs(s).plph(myf(1:5:end),:,chi,i))));
%                     v = imag(exp(sqrt(-1)*sq(OutArgs(s).plph(myf(1:5:end),:,chi,i))));
%                     quiver(phbin,f(myf(1:5:end)),u,v,1, 'Color','k');


                    %         subplot(ni,nj, 6)
                    %         imagesc(phbin,f(myf),abs(sq(OutArgs(s).ycom(:,myf,1,2)))');axis xy
                    %         title('average coherence spks to eeg');
                    %         hold on
                    %
                    %         subplot(ni,nj, 7)
                    %         imagesc(phbin,f(myf),(sq(OutArgs(s).phinorm(:,myf,1,2)))');axis xy
                    %         title('norm phase of spks to eeg');
                    %
                    %         subplot(ni,nj, 8)
                    %         imagesc(phbin,f(myf),(sq(OutArgs(s).phi(:,myf,1,2)))');axis xy
                    %         title('phase of spks to eeg');
                    %

                    tit = [FileBase ' : ' num2str(OutArgs(s).ElClu(i,:)) ...
                        ' (' OutArgs(s).State '), Rate=' num2str(OutArgs(s).FirRate(i))];
                    %h = suptitle(tit);
                    %set(h,'FontSize',4);
                    subplot(ni,nj,1);
                    h = title(tit);
                    set(h,'FontSize',4);
                    [x y b] = PointInput(1);
                    if b==2
                        return;
                    elseif b==3
                        chi = round(y);
                    elseif b==1
                        break;
                    end
                end % of while
            end

        end
end





