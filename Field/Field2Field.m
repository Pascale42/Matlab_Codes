function out = Field2Field(FileBase,varargin)
%function Field2Field(FileBase,fMode, FreqRange,State,Window)
[fMode, FreqRange,States,Window] = DefaultArgs(varargin,{'compute',[1 150], {'REM','SWS','RUN'},2^10});

Par = LoadPar([FileBase '.xml']);
MinPeriod = 3; %seconds for theta periods selection
if ischar(States); States = {States}; end;
nSt = length(States);

[RepCh ifTetr Info] = RepresentChan(FileBase);
nRepCh = length(RepCh);
eFs = 1250;
switch fMode
    case  'compute'
        for w=1:nSt
            Period = SelectStates(FileBase, States{w},eFs*2);
            if isempty(Period)
                fprintf('No %s periods. Return empty output\n',States{w});
                continue;
            end
            eeg = LoadBinary([FileBase '.eeg'],RepCh,Par.nChannels,Period);
            eeg = WhitenSignal(eeg);
            [out(w).y out(w).f out(w).phi out(w).yerr, out(w).phierr, out(w).phloc] = ...
                mtptchd(eeg,[],[],2^nextpow2(Window),eFs,2^nextpow2(Window),[],2,'linear',[],FreqRange);
            out(w).State = States{w};
            out(w).Channels = RepCh;
        end

    case 'display'
        load([FileBase '.' mfilename '.mat']);
        for ii=1:length(OutArgs)
            if isempty(OutArgs(ii).State)
                OutArgs(ii)=[];
            end
        end

        nStates = length(OutArgs);

        chi=1; fi =round(length(OutArgs(1).f)/2);
        nCh = length(OutArgs(1).Channels);

        nx=3; ny=3;
        clickmode=1;
        redraw = 1;
        curplot=3;
        % axlayout = reshape(
        while 1
            if redraw
                fprintf('chi = %d, fi = %d, curplot = %d',chi,fi,curplot);
                myi = setdiff([1:nCh],chi);
                chnum = OutArgs(1).Channels(chi);
                for s=1:nStates
                    
                    ych = rem(chnum,Info.nChInShank);
                    xch = ceil((chnum+eps)/Info.nChInShank);
                    mycoh = []; myph=[]; mypow=[];
                    mycoh(:,myi) = atanh(sq(OutArgs(s).y(:,chi,myi)));
                    mycoh(:,chi) = NaN*ones(size(mycoh,1),1);
                    myph(:,myi) = (sq(OutArgs(s).phi(:,chi,myi)));
                    myph(:,chi) = NaN*ones(size(mycoh,1),1);

                    for jj=1:nCh
                        mypow(:,jj) = 10*log10(sq(OutArgs(s).y(:,jj,jj)));
                    end
                    figure(100+s);
                    vars = {'mypow','mycoh','myph'};
                    lab ={'power','coherence','phase'};
                    for k=1:3
                        mat = eval(vars{k});
                        subplot(ny,nx,(k-1)*nx+[1 2])
                        imagesc(OutArgs(s).f,[1:nCh],mat');hold on
                        % colorbar('NorthOutside');
                        title([lab{k} ' : ' OutArgs(s).State]);
                        ylabel('channels'); xlabel('Freq');
                        hold on
                        plot(OutArgs(s).f(fi),chi, 'ko','MarkerSize',10);
                        Lines(OutArgs(s).f(fi),[], 'k');
                        hold off
                        if k==3
                            caxis([-pi pi]);
                        end

                        
                        subplot(ny,nx,(k-1)*nx+[3])
                        MapSilicon(mat(fi,:),OutArgs(s).Channels);
                        hold on
                        plot(xch, ych, 'ko','MarkerSize',10);
                        hold off
                        if k==3
                            caxis([-pi pi]);
                        end

                    end
                end %end state loop
            end %if redraw
            figure(100+1);
            chld = get(gcf,'Children');
            chld = flipud(chld(:));
            set(gcf,'CurrentAxes',chld(curplot));
%             if exist('axH','var')
%                 set(0,'CurrentFigure',figH);
%                 set(figH,'CurrentAxes',axH);
%             end
            [x y b key] = PointInput(1);
            if x<0 | y<0
                fprintf('bad point\n');
                [x y b] = PointInput(1);
            end
           % x = x(2); y=y(2); b=b(2);
%           [myax myfig axH figH] = WhereIClicked;
            
            if b==1 | b==3
                               
                if rem(curplot,2)>0%rem(myax,2)>0 %so it is 1 or 2 -chan x freq plot
%                       clickmode=1;

                       [dummy fi] = min(abs(OutArgs(s).f-x));
                    if b==1                       
                       chi = round(y); chi = max(chi,1); chi = min(chi,size(OutArgs(s).y,2));
                    end
                 
                else
                   % clickmode=2;
                    x=round(x);
                    x= min(max(x,1),Info.nShanks);
                    y = round(y);
                    y= min(max(y,1),Info.nChInShank);
                    chnum = (x-1)*Info.nChInShank+y;
                    [dummy chi] = min(abs(OutArgs(s).Channels-chnum));
                end
                redraw = 1;
            elseif b==0
                
                %which column to focus on
                if double(key)>=49 & double(key)<60
                    curplot = str2num(key);
                end
                redraw = 0;
%                 if clickmode==1
%                     switch double(key)
%                         case 28
%                             fi=fi-1; fi=max(fi,1);
%                         case 29
%                             fi=fi+1; fi=min(fi,max(OutArgs(s).f));
%                         case 30
%                             chi =chi-1; chi = max(1,chi);
%                         case 31
%                             chi =chi+1; chi = min(chi,nCh);
%                     end
%                 else
%                      switch double(key)
%                         case 28
%                              xch=xch-1; xch = max(1,xch);
%                          case 29
%                              xch=xch+1; xch = min(xch,Info.nShanks);
%                          case 30
%                               ych=ych-1; ych = max(1,ych);
%                          case 31
%                              ych=ych+1; ych = min(ych,Info.nChInShank);
%                              
%                      end
% 
%                     chnum = (x-1)*Info.nChInShank+y;
%                     [dummy chi] = min(abs(OutArgs(1).Channels-chnum));  
%                 end

            else
                return;
            end

        end %end while

end