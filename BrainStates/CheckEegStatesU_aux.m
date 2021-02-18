function ChechEegStatesU_aux(action)

global gCheckEegStates

%try

    switch action

        case 'keyboard'

            %check what key was pressed in figure
            whatkey = get(gcf,'CurrentCharacter');
            %fprintf('key pressed: %d\n',double(whatkey));
               switch double(whatkey)
                    
                    case double('h') %help

                        msgbox({'The following keys switch (on/off) the following behavior:', ...
                            't (default) - to move the current position pointer and display the corresponding trace of eeg', ...
                            'n - add new region first click set''s one border, second - the other. new region is displayed as green/red lines',...
                            'm - move one border. click next to the border and drag it where you want. ',...
                            'd - delete the border. Just click on the border',...
                            'z - toggles zoom state. Mouse :left/right btn - zoom in/out the position of the mouse. Keybd: f - resets x axis to max. d/i - zoom in/out the current position',...
                            'f - and then  up/down - zoom in/out the y axis in spectrograms',...
                            'w - and then up/down - decrease/increas the size of the window for the traces display ',...
                            'u - update periods. Removes deleted and reorders all lines, pairs consequituive ones. Colorizes blue -beg, red -end',...
                            's - save the results in the file (currently ASCII two columnt, 1250 s.r.', ...
                            'arrows left/right - move in the spectrogram' });

                    case double('t')
                        keyboard
                        set(gcf,'Pointer','arrow');
                        set(gCheckEegStates.figh, 'Name', 'CheckEegStates : traces');
                        gCheckEegStates.Mode = 't';
                        set(gcf,'WindowButtonMotionFcn','', 'WindowButtonDown','CheckEegStatesU_aux(''mouseclick'')','WindowButtonUpFcn','');

                    case 32 %space
                        set(gcf,'Pointer','arrow');
                        set(gCheckEegStates.figh, 'Name', 'CheckEegStates : traces');
                        gCheckEegStates.Mode = 't';
                        set(gcf,'WindowButtonMotionFcn','', 'WindowButtonDown','CheckEegStatesU_aux(''mouseclick'')','WindowButtonUpFcn','');

                    case double('z') % zoom in x axis
                        % stupid matlab7 zoom redefines all my WindowButtonXXX functions  -
                        % so have to go back to matlab6.5 zoom function - hacked as
                        % oldzoom in my matlab/draft. maybe some day I will write my
                        % own - but for now it works fine. wrote my own!!! no more
                        % zoom!!
                        MyPointer('eye');
                        %oldzoom xon
                        gCheckEegStates.Mode='z';
                        set(gCheckEegStates.figh, 'Name', 'CheckEegStates : zoom');
                        %set(gcf,'WindowButtonUpFcn','CheckEegStatesU_aux(''zoomequal'')');
                        set(gcf,'KeyPressFcn','CheckEegStatesU_aux(''zoomkeys'')');

                    case double('n') % add new region
                        if strcmp(gCheckEegStates.Mode, 'n')
                            gCheckEegStates.Mode = 't';
                            set(gcf,'Pointer','arrow');
                            set(gCheckEegStates.figh, 'Name', 'CheckEegStates : traces');
                            %set(gcf,'WindowButtonMotionFcn','');
                        else
                            MyPointer('pencil');
                            set(gCheckEegStates.figh, 'Name', 'CheckEegStates : add new region');
                            gCheckEegStates.Mode = 'n';
                            %set(gCheckEegStates.figh,'WindowButtonMotionFcn','CheckEegStatesU_aux(''traces'')');
                        end


                    case double('q') % quit
                        close(gCheckEegStates.figh);
                        clear global gCheckEegStates;
                        return
                        %set(gcf,'WindowButtonDownFcn','');
                        %set(gcf,'KeyPressFcn', '');

                    case double('s') % save periods
                        CheckEegStatesU_aux('update_per');
                        if ~isempty(gCheckEegStates.Periods)
                            SaveFileName = [gCheckEegStates.FileBase '.sts.' gCheckEegStates.State];
                            saveans = inputdlg({'Enter filename '},'Save Periods',1,{SaveFileName});
                            if ~strcmp(saveans,'No')
                                msave(saveans{1}, round(gCheckEegStates.Periods*1250));

                                %MaveEvtFiles(CheckEegStates.Periods, [saveans{1} '.evt']);
                            end
                        end

                    case 29 % move forward - right arrow
             
                        gCheckEegStates.Mode = 't';
                        child = get(gcf,'Children');
                        curxlim = get(child(end),'XLim');
                        step = diff(curxlim)/4;
                        gCheckEegStates.t = min(gCheckEegStates.trange(2), gCheckEegStates.t + step);
                        for ii=1:gCheckEegStates.nPlots-1
                            subplot(gCheckEegStates.nPlots,1,ii);
                            newax(1) = max(gCheckEegStates.trange(1), curxlim(1)+step);
                            newax(2) = min(gCheckEegStates.trange(2), curxlim(2)+step);
                            xlim(newax);
                        end
                        CheckEegStatesU_aux('traces');


                    case 28  % move backward - left arrow
                        gCheckEegStates.Mode = 't';
                        child = get(gcf,'Children');
                        curxlim = get(child(end),'XLim');
                        step = diff(curxlim)/4;
                        gCheckEegStates.t = max(gCheckEegStates.trange(1), gCheckEegStates.t - step);
                        for ii=1:gCheckEegStates.nPlots-1
                            subplot(gCheckEegStates.nPlots,1,ii);
                            newax(1) = max(gCheckEegStates.trange(1), curxlim(1)-step);
                            newax(2) = min(gCheckEegStates.trange(2), curxlim(2)-step);
                            xlim(newax);
                        end
                        CheckEegStatesU_aux('traces');


                    case double('w') % change the window size
                        gCheckEegStates.Mode = 'w';
                        set(gCheckEegStates.figh, 'Name', 'CheckEegStates : trace window resize');
                        set(gCheckEegStates.figh,'KeyPressFcn', 'CheckEegStatesU_aux(''windowsize'')');
                        %CheckEegStatesU_aux('windowsize');

                    case double('f') % change the freq. range
                        gCheckEegStates.Mode = 'f';
                        set(gCheckEegStates.figh, 'Name', 'CheckEegStates : freq. axis resize');
                        set(gCheckEegStates.figh,'KeyPressFcn', 'CheckEegStatesU_aux(''freqsize'')');

                    case double('u') % update periods borders coloring etc

                        CheckEegStatesU_aux('update_per');

                    case 'c' %color%cycle through the periods
                        gCheckEegStates.Mode = 'c';
                        set(gCheckEegStates.figh, 'Name', 'CheckEegStates : color adjust');
                        set(gCheckEegStates.figh,'KeyPressFcn', 'CheckEegStatesU_aux(''coloradj'')');
                        %to be done

                end
                if ~isempty(gCheckEegStates.Periods)
                    switch whatkey

                        case 'm' % move regions
                            if strcmp(gCheckEegStates.Mode, 'm')
                                gCheckEegStates.Mode = 't';
                                set(gcf,'Pointer','arrow');
                                set(gCheckEegStates.figh, 'Name', 'CheckEegStates : traces');
                            else
                                MyPointer('mover');
                                set(gCheckEegStates.figh, 'Name', 'CheckEegStates : move line');
                                gCheckEegStates.Mode = 'm';
                            end


                        case 'd' % delete regions
                            if strcmp(gCheckEegStates.Mode, 'd')
                                gCheckEegStates.Mode = 't';
                                set(gcf,'Pointer','arrow');
                                set(gCheckEegStates.figh, 'Name', 'CheckEegStates : traces');
                            else
                                MyPointer('scull');
                                set(gCheckEegStates.figh, 'Name', 'CheckEegStates : delete line');
                                gCheckEegStates.Mode = 'd';
                            end

                        case 'e' %erase all periods

                            gCheckEegStates.Periods = [];
                            for ii=1:gCheckEegStates.nPlots
                                subplot(gCheckEegStates.nPlots,1,ii);
                                notnan = find(~isnan(gCheckEegStates.lh{ii}(:)) & gCheckEegStates.lh{ii}(:)>0 );
                                delete(gCheckEegStates.lh{ii}(notnan));
                                gCheckEegStates.lh{ii}=[];
                            end


                    end % for switch whatkey
                end % for if


        case 'mouseclick'
            
            whatbutton = get(gcf,'SelectionType');
            mousecoord = get(gca,'CurrentPoint');
            xmouse = mousecoord(1,1);
            curaxis = get(gcf,'CurrentAxes');
 %matlab has stupid doubl click hadnling - doesnt work!!!           
% %            now need to catch double click
%             if strcmp(whatbutton,'open') 
%                 buttontype = gCheckEegStates.LastBut;
%             else
%                 buttontype = whatbutton;
%             end
%             gCheckEegStates.LastBut = whatbutton;

%            fprintf('%s %s\n',whatbutton,buttontype);
            %switch gCheckEegStates.Mode
            %   case 't' % default
            if strcmp(gCheckEegStates.Mode,'d')

                if sum(~isnan(gCheckEegStates.Periods(:)))==0
                    gCheckEegStates.Mode = 't';
                    set(gcf,'Pointer','arrow');
                    set(gCheckEegStates.figh, 'Name', 'CheckEegStates : traces');
                else
                    %find closest line to clicked point
                    [dist2ln(1) closest_ln(1)] = min(abs(gCheckEegStates.Periods(:,1)-xmouse));
                    [dist2ln(2) closest_ln(2)] = min(abs(gCheckEegStates.Periods(:,2)-xmouse));
                    [dummy whichone] = min(dist2ln); % choses which one is closer left or right
                    closest_ln = closest_ln(whichone); %that gives the period index
                    gCheckEegStates.SelLine = [closest_ln whichone]; %store for outside use
                    for ii=1:gCheckEegStates.nPlots
                        subplot(gCheckEegStates.nPlots,1,ii);
                        delete(gCheckEegStates.lh{ii}(closest_ln,whichone)); %delete particular line
                        gCheckEegStates.lh{ii}(closest_ln,whichone) = NaN; %set the handle to undefined
                    end
                    gCheckEegStates.Periods(closest_ln,whichone)=NaN;

                end

                %case 'z' % my own zoom with mouse - damn matlab zoom!
            %    elseif strcmp(gCheckEegStates.Mode,'z')
            elseif (strcmp(gCheckEegStates.Mode,'z')&& (strcmp(whatbutton,'normal') || strcmp(whatbutton,'alt'))) || strcmp(whatbutton,'open')
            %elseif strcmp(whatbutton,'open')
                whatbutton = get(gcf,'SelectionType');
                mousecoord = get(gca,'CurrentPoint');
                xmouse = mousecoord(1,1);

                switch whatbutton

                    case 'normal'
                        factor = 1/2;
                    case 'alt'
                        factor= 2;
                    otherwise
                        factor= 1/2;
                        %do nothing

                end
                if xmouse>gCheckEegStates.trange(1) & xmouse < gCheckEegStates.trange(2)
                    gCheckEegStates.t = xmouse; % set current poisiton to the clicked point
                end
                if exist('factor','var')
                    for ii=1:gCheckEegStates.nPlots-1
                        subplot(gCheckEegStates.nPlots,1,ii);
                        curax = get(gca,'XLim');
                        curcenter =gCheckEegStates.t;
                        curwidth = diff(curax);
                        newax(1) = max(gCheckEegStates.trange(1), curcenter-curwidth/2*factor);
                        newax(2) = min(gCheckEegStates.trange(2), curcenter+curwidth/2*factor);
                        xlim(newax);
                    end
                end
                CheckEegStatesU_aux('traces');

                %case 'n' % add new region

            elseif (strcmp(gCheckEegStates.Mode,'n')&&strcmp(whatbutton,'normal')) || strcmp(whatbutton,'alt') % right click or click+ctrl
                %gCheckEegStates.Mode = 'n';
                MyPointer('pencil');
                set(gCheckEegStates.figh, 'Name', 'CheckEegStates : add new region');
                %set(gCheckEegStates.figh,'WindowButtonMotionFcn','CheckEegStatesU_aux(''traces'')');
                %                    switch whatbutton
                %case 'normal'
                gCheckEegStates.Periods = [gCheckEegStates.Periods; [xmouse NaN]];% add new period
                fprintf('Added period at %2.2f seconds\n',xmouse);
                lnindex = size(gCheckEegStates.Periods,1);
                %plot new lines
                for ii=1:gCheckEegStates.nPlots
                    subplot(gCheckEegStates.nPlots,1,ii);
                    gCheckEegStates.lh{ii}(lnindex,1) = Lines(xmouse,[],'b',[],2);
                end

                %                         case 'alt'
                %
                %                             if isempty(gCheckEegStates.newl)
                %                                 gCheckEegStates.newl(1) = xmouse; %first border
                %
                %                             else
                %                                 gCheckEegStates.newl(2) = xmouse; % second border
                %                                 gCheckEegStates.newl = sort(gCheckEegStates.newl); %order
                %                                 gCheckEegStates.Periods = [gCheckEegStates.Periods; gCheckEegStates.newl(:)'];% add new period
                %                                 fprintf('Added period %2.2f to %2.2f seconds\n',gCheckEegStates.newl(1),gCheckEegStates.newl(2));
                %                                 lnindex = size(gCheckEegStates.Periods,1);
                %                                 %plot new lines
                %                                 for ii=1:gCheckEegStates.nPlots
                %                                     subplot(gCheckEegStates.nPlots,1,ii);
                %                                     gCheckEegStates.lh{ii}(lnindex,1) = Lines(gCheckEegStates.newl(1),[],'b',[],2);
                %                                     gCheckEegStates.lh{ii}(lnindex,1) = Lines(gCheckEegStates.newl(2),[],'r',[],2);
                %                                 end
                %                                 gCheckEegStates.newl = [];
                %                                 %CheckEegStatesU_aux('lines'); %plot new lines in all plots
                %
                %                             end

                %    end


                %case 'm' % move lines
            elseif (strcmp(gCheckEegStates.Mode,'m') &&strcmp(whatbutton,'normal'))|| strcmp(whatbutton,'extend') % middle click or click+shift
                MyPointer('mover');
                set(gCheckEegStates.figh, 'Name', 'CheckEegStates : move line');
                %gCheckEegStates.Mode = 'm';

                %find closest line to clicked point
                [dist2line(1) closest_ln(1)] = min(abs(gCheckEegStates.Periods(:,1)-xmouse));
                [dist2line(2) closest_ln(2)] = min(abs(gCheckEegStates.Periods(:,2)-xmouse));
                [dummy whichone] = min(dist2line); % choses which one is closer: left or right
                closest_ln = closest_ln(whichone); %that gives the period index
                gCheckEegStates.SelLine = [closest_ln whichone]; %store for outside use
                fprintf('Moved %d line of period # %d from %2.2f',whichone,closest_ln,gCheckEegStates.Periods(closest_ln,whichone));
                curax = get(gcf,'CurrentAxes');

                for ii=1:gCheckEegStates.nPlots
                    subplot(gCheckEegStates.nPlots,1,ii);
                    set(gCheckEegStates.lh{ii}(closest_ln,whichone),'LineWidth',2);
                end
                set(gcf,'CurrentAxes',curax);
                %setup callbacks
                set(gcf,'WindowButtonMotionFcn','CheckEegStatesU_aux(''linemove'')'); % move the line along the
                set(gcf,'WindowButtonUpFcn','CheckEegStatesU_aux(''linemove_end'')'); %once relieased - change the coordinates of the lines and update plots

                %case 'd' %delete

            elseif (strcmp(gCheckEegStates.Mode,'t') && strcmp(whatbutton,'normal')) || strcmp(whatbutton,'normal')
                set(gcf,'WindowButtonMotionFcn','', 'WindowButtonDown','CheckEegStatesU_aux(''mouseclick'')','WindowButtonUpFcn','');
                gCheckEegStates.Mode = 't';
                set(gcf,'Pointer','arrow');
                set(gCheckEegStates.figh, 'Name', 'CheckEegStates : traces');
                set(gcf,'WindowButtonMotionFcn','');

                gCheckEegStates.t = xmouse;

                % check if the browising point got out of the spectrograms zoom
                % window  - can happen if click on the traces window. adjust
                subplot(gCheckEegStates.nPlots,1,1);
                curxlim = xlim;
                curwidth = diff(xlim);
                if gCheckEegStates.t > curxlim(2)-gCheckEegStates.Window/1250%0.1*curwidth
                    shift = diff(curxlim)-gCheckEegStates.Window/1250;
                elseif gCheckEegStates.t < curxlim(1)+gCheckEegStates.Window/1250%0.1*curwidth
                    shift = -diff(curxlim)+gCheckEegStates.Window/1250;
                end
                if exist('shift','var')
                    for ii=1:gCheckEegStates.nPlots-1
                        subplot(gCheckEegStates.nPlots,1,ii);
                        newax(1) = max(gCheckEegStates.trange(1), curxlim(1)+shift);
                        newax(2) = min(gCheckEegStates.trange(2), curxlim(2)+shift);
                        xlim(newax);
                    end
                end
                CheckEegStatesU_aux('traces');


            end

        case 'traces'

            if 0 % this is for online trace update in the mode of adding new periods
                if strcmp(gCheckEegStates.Mode,'n')
                    mousecoord = get(gca,'CurrentPoint');
                    xmouse = mousecoord(1,1);
                    gCheckEegStates.t = xmouse;
                end
            end

            %plot current position of the cursor

            for ii=1:gCheckEegStates.nPlots-1
                subplot(gCheckEegStates.nPlots,1,ii);
                if ~isempty(gCheckEegStates.cposh{ii})
                    delete(gCheckEegStates.cposh{ii});
                end
                gCheckEegStates.cposh{ii} = Lines(gCheckEegStates.t,[],'k','--',2);
            end

            %load the traces from  eeg file
            maxlen = FileLength([gCheckEegStates.FileBase '.eeg'])/2/gCheckEegStates.nChannels;
            segbeg = max(1,round(gCheckEegStates.t*1250-gCheckEegStates.Window*3/2));
            seglen = min(round(gCheckEegStates.Window*3), maxlen-segbeg-1);
            Seg = LoadSegs([gCheckEegStates.FileBase '.eeg'],segbeg,seglen,gCheckEegStates.nChannels,gCheckEegStates.Channels,1250);

            % plot them in the lowest subplot
            figure(gCheckEegStates.figh);
            subplot(gCheckEegStates.nPlots,1,gCheckEegStates.nPlots)
            cla
            htr = PlotTraces(Seg,(segbeg+[0:seglen-1])/1250,1250,3);
            ylabel(num2str(fliplr(gCheckEegStates.Channels(:)')));
            set(gca,'YTick',[]);

            %plot the lines in traces subplot
            if ~isempty(gCheckEegStates.Periods)
                gCheckEegStates.lh{gCheckEegStates.nPlots}(:,1) = Lines(gCheckEegStates.Periods(:,1),[],'b',[],2);
                gCheckEegStates.lh{gCheckEegStates.nPlots}(:,2) = Lines(gCheckEegStates.Periods(:,2),[],'r',[],2);
            end
            gCheckEegStates.cposh{gCheckEegStates.nPlots} = Lines(gCheckEegStates.t+[-1 1]*gCheckEegStates.tstep,[],'k','--',2);

            % now plot cool lines to indicqate where the traces come from in
%             % spectrograms - looks great, ah? :)
%             if ~isempty(gCheckEegStates.coolln)
%                 delete(gCheckEegStates.coolln);
%             end
%             %get the coordinates
%             subplot(gCheckEegStates.nPlots,1,gCheckEegStates.nPlots-2);
%             ax0 = get(gca,'Position');
%             xlim0 = get(gca,'XLim');
%             
%             subplot(gCheckEegStates.nPlots,1,gCheckEegStates.nPlots-1);
%             ax1 = get(gca,'Position');
%             xlim1 = get(gca,'XLim');
%             
%             subplot(gCheckEegStates.nPlots,1,gCheckEegStates.nPlots);
%             ax2 = get(gca,'Position');
%             xlim2 = get(gca,'XLim');
%             
%             lnx = [repmat(ax0(1)+ax0(3)*(gCheckEegStates.t-xlim0(1))/diff(xlim0),1,2); ...
%                 repmat(ax2(1)+ax2(3)*(gCheckEegStates.t-xlim2(1))/diff(xlim2),1,2)+ax2(3)*[-1 1]*gCheckEegStates.tstep/diff(xlim2)];
%             lny = repmat([ax1(2)  ax2(2)+ax2(4)]',1,2);
%             %plot the lines
%             gCheckEegStates.coolln(1) = annotation(gcf, 'line', lnx(:,1), lny(:,1));
%             gCheckEegStates.coolln(2) = annotation(gcf, 'line', lnx(:,2), lny(:,2));


        case 'lines'
            if ~isempty(gCheckEegStates.Periods)
                % plot new lines
                for ii=1:gCheckEegStates.nPlots
                    subplot(gCheckEegStates.nPlots,1,ii);
                    if ~isempty(gCheckEegStates.lh{ii}(:))
                        gi = find(~isnan(gCheckEegStates.lh{ii}(:)));
                        delete(gCheckEegStates.lh{ii}(gi)); %delete all existing lines - if not deleted already
                        gCheckEegStates.lh{ii}=[];
                    end
                    gi = find(~isnan(gCheckEegStates.Periods(:,1)));
                    gCheckEegStates.lh{ii}(gi,1) = Lines(gCheckEegStates.Periods(gi,1),[],'b',[],2);
                    gi = find(~isnan(gCheckEegStates.Periods(:,2)));
                    gCheckEegStates.lh{ii}(gi,2) = Lines(gCheckEegStates.Periods(gi,2),[],'r',[],2);
                end
            end

        case 'linemove'
            whichln = gCheckEegStates.SelLine;
            mousecoord = get(gca,'CurrentPoint');
            xmouse = mousecoord(1,1);
            col = {'b','r'};
            curax = get(gcf,'CurrentAxes');
            for ii=1:gCheckEegStates.nPlots
                subplot(gCheckEegStates.nPlots,1,ii);
                set(gCheckEegStates.lh{ii}(whichln(1),whichln(2)),'XData',[1 1]*xmouse);
            end
            set(gcf,'CurrentAxes',curax);


        case 'linemove_end' % callback for the butto-nup of the line move action

            mousecoord = get(gca,'CurrentPoint');
            xmouse = mousecoord(1,1);
            whichln = gCheckEegStates.SelLine;
            col = {'b','r'};
            for ii=1:gCheckEegStates.nPlots
                subplot(gCheckEegStates.nPlots,1,ii);
                %delete(gCheckEegStates.lh{ii}(whichln(1), whichln(2))); %delete particular line
                %gCheckEegStates.lh{ii}(whichln(1),whichln(2)) = Lines(xmouse,[],col{whichln(2)});
                set(gCheckEegStates.lh{ii}(whichln(1),whichln(2)),'XData',[1 1]*xmouse);
            end
            gCheckEegStates.Periods(whichln(1),whichln(2)) = xmouse; %update the periods matrix
            fprintf(' to %2.2f\n',gCheckEegStates.Periods(whichln(1),whichln(2)));
            set(gcf,'WindowButtonUpFcn','');
            set(gcf,'WindowButtonMotionFcn','');


            % case 'zoomequal' % makes all axes equal when zooming in one
            % 	%childax = get(gcf,'Children')
            %     %gca
            % 	%if ~ismember(get(gcf,'CurrentAxes'),child)
            % 	%	set(gcf,'CurrentAxes',child(1));
            % 	%end
            %     %drawnow
            % 	curaxis = get(gca,'XLim')
            % 	for ii=1:gCheckEegStates.nPlots-1
            %     	subplot(gCheckEegStates.nPlots,1,ii);
            % 		xlim(curaxis(1:2));
            %     end
            % 	% also update pointer and traces if outside the range
            %     gCheckEegStates.t
            %     if gCheckEegStates.t > curaxis(2) | gCheckEegStates.t < curaxis(1)
            %            gCheckEegStates.t = mean(curaxis);
            %     end
            %     CheckEegStatesU_aux('traces');

        case 'windowsize' %kbd callback for the traces window size
            whatkey = get(gcf,'CurrentCharacter');
            switch double(whatkey)
                case 30
                    gCheckEegStates.Window = gCheckEegStates.Window*2;
                    CheckEegStatesU_aux('traces');
                case 31
                    gCheckEegStates.Window = gCheckEegStates.Window/2;
                    CheckEegStatesU_aux('traces');
                case double('w') %return to base mode
                    gCheckEegStates.Mode = 't';
                    set(gCheckEegStates.figh, 'Name', 'CheckEegStates : traces');
                    set(gCheckEegStates.figh,'KeyPressFcn', 'CheckEegStatesU_aux(''keyboard'')');

                otherwise
                    fprintf('wrong key: use i to increase\n d to decrease\n w - to return to display mode');

            end


        case 'freqsize' % kbd functionality to change freq. axes

            whatkey = get(gcf,'CurrentCharacter');
            switch double(whatkey)
                case  30
                    factor = 2;

                case 31
                    factor = 1/2;

                case double('f')
                    gCheckEegStates.Mode = 't';
                    set(gCheckEegStates.figh, 'Name', 'CheckEegStates : traces');
                    set(gCheckEegStates.figh,'KeyPressFcn', 'CheckEegStatesU_aux(''keyboard'')');

                otherwise
                    fprintf('wrong key: use i to increase\n d to decrease\n w - to return to display mode\n');

            end

            if exist('factor','var')
                for ii=1:gCheckEegStates.nPlots-2
                    subplot(gCheckEegStates.nPlots,1,ii);
                    if ii==1
                        curax = get(gca,'YLim');
                        newax = [1 min(gCheckEegStates.FreqRange(2), curax(2)*factor)];
                    end
                    ylim(newax);
                end
            end
        
        case 'coloradj' % kbd functionality to color adjust

            whatkey = get(gcf,'CurrentCharacter');
            switch double(whatkey)
                case  30
                    factor = -0.1;

                case 31
                    factor = +0.1;

                case double('c')
                    gCheckEegStates.Mode = 't';
                    set(gCheckEegStates.figh, 'Name', 'CheckEegStates : traces');
                    set(gCheckEegStates.figh,'KeyPressFcn', 'CheckEegStatesU_aux(''keyboard'')');

                otherwise
                    fprintf('wrong key: use i to increase\n d to decrease\n w - to return to display mode\n');

            end

            if exist('factor','var')
                for ii=1:gCheckEegStates.nPlots-2
                    subplot(gCheckEegStates.nPlots,1,ii);
                    if ii==1
                        curax = get(gca,'CLim');
                        newax = [curax(1) curax(2)*(1+factor)];
                    end
                    caxis(newax);
                end
            end

            
        case  'zoomkeys' % kbd callback for the zoom
            whatkey = get(gcf,'CurrentCharacter');
            switch double(whatkey)
                case  30 % up
                    factor = 2;

                case 31 %down
                    factor = 1/2;

                case double('f')
                    for ii=1:gCheckEegStates.nPlots-1
                        subplot(gCheckEegStates.nPlots,1,ii);
                        xlim(gCheckEegStates.trange);
                    end

                case double('z') %return to base mode
                    oldzoom off
                    set(gCheckEegStates.figh, 'Name', 'CheckEegStates : traces');
                    set(gcf,'WindowButtonUpFcn','');
                    set(gcf,'KeyPressFcn','CheckEegStatesU_aux(''keyboard'')');
                    gCheckEegStates.Mode = 't';
                    set(gcf,'Pointer','arrow');

                otherwise
                    fprintf('wrong key: use i to increase\nd to decrease\nz - to return to display mode\n');

            end

            if exist('factor','var')
                for ii=1:gCheckEegStates.nPlots-1
                    subplot(gCheckEegStates.nPlots,1,ii);
                    curax = get(gca,'XLim');
                    %  curcenter = mean(curax);
                    curcenter =gCheckEegStates.t;
                    curwidth = diff(curax);
                    newax(1) = max(gCheckEegStates.trange(1), curcenter-curwidth/2*factor);
                    newax(2) = min(gCheckEegStates.trange(2), curcenter+curwidth/2*factor);
                    xlim(newax);
                end
            end

        case 'update_per' % sort out the periods (eliminate NaN and reorder), colorize according to left/right

            %fisrt check length
            if size(gCheckEegStates.Periods,1)~=size(gCheckEegStates.lh{1},1)
                error('periods array and line handles arrays are different in length!');
            end


            %strech and remove NaN

            newPeriods = gCheckEegStates.Periods(:);
            notnan = find(~isnan(newPeriods));
            if ~mod(length(notnan),2)
                [sortedPeriods sortind] = sort(newPeriods(notnan));
                gCheckEegStates.Periods = reshape(sortedPeriods,2,[])';
                %update the handles
                for ii=1:gCheckEegStates.nPlots
                    subplot(gCheckEegStates.nPlots,1,ii);
                    newlh = gCheckEegStates.lh{ii}(:);
                    newlh = newlh(notnan);
                    gCheckEegStates.lh{ii} = reshape(newlh(sortind),2,[])';
                    set(gCheckEegStates.lh{ii}(:,1),'Color','b');
                    set(gCheckEegStates.lh{ii}(:,2),'Color','r');
                end
            else
                fprintf('one border is missing, check and update again\n');
            end


    end
% catch
%     set(gcf,'WindowButtonUpFcn','');
%     set(gcf,'WindowButtonMotionFcn','');
%     set(gcf,'WindowButtonDownFcn','');
%     lasterr
%     keyboard
% 
% end