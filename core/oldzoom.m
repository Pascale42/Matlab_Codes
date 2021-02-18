function out = oldzoom(varargin)
%OLDZOOM   Zoom in and out on a 2-D plot.
%   OLDZOOM with no arguments toggles the oldzoom state.
%   OLDZOOM(FACTOR) oldzooms the current axis by FACTOR.
%       Note that this does not affect the oldzoom state.
%   OLDZOOM ON turns oldzoom on for the current figure.
%   OLDZOOM XON or OLDZOOM YON turns oldzoom on for the x or y axis only.
%   OLDZOOM OFF turns oldzoom off in the current figure.
%   OLDZOOM RESET resets the oldzoom out point to the current oldzoom.
%   OLDZOOM OUT returns the plot to its current oldzoom out point.
%   If OLDZOOM RESET has not been called this is the original
%   non-oldzoomed plot.  Otherwise it is the oldzoom out point
%   set by OLDZOOM RESET.
%
%   When oldzoom is on, click the left mouse button to oldzoom in on the
%   point under the mouse.  Click the right mouse button to oldzoom out.
%   Each time you click, the axes limits will be changed by a factor 
%   of 2 (in or out).  You can also click and drag to oldzoom into an area.
%   It is not possible to oldzoom out beyond the plots' current oldzoom out
%   point.  If OLDZOOM RESET has not been called the oldzoom out point is the
%   original non-oldzoomed plot.  If OLDZOOM RESET has been called the oldzoom out
%   point is the oldzoom point that existed when it was called.
%   Double clicking oldzooms out to the current oldzoom out point - 
%   the point at which oldzoom was first turned on for this figure 
%   (or to the point to which the oldzoom out point was set by OLDZOOM RESET).
%   Note that turning oldzoom on, then off does not reset the oldzoom out point.
%   This may be done explicitly with OLDZOOM RESET.
%   
%   OLDZOOM(FIG,OPTION) applies the oldzoom command to the figure specified
%   by FIG. OPTION can be any of the above arguments.

%   OLDZOOM FILL scales a plot such that it is as big as possible
%   within the axis position rectangle for any azimuth and elevation.

%   Clay M. Thompson 1-25-93
%   Revised 11 Jan 94 by Steven L. Eddins
%   Copyright 1984-2002 The MathWorks, Inc. 
%   $Revision: 5.63 $  $Date: 2002/04/08 21:44:39 $

%   Note: oldzoom uses the figure buttondown and buttonmotion functions
%
%   OLDZOOM XON oldzooms x-axis only
%   OLDZOOM YON oldzooms y-axis only


%
% PARSE ARGS - set fig and oldzoomCommand
%
switch nargin
    % no arg in
case 0
    oldzoomCommand='toggle';
    fig=get(0,'currentfigure');
    if isempty(fig)
        return
    end
    % one arg in 
case 1
    % If argument is string, it is a oldzoom command
    % (i.e. (on, off, down, xdown, etc.). 
    if isstr(varargin{1})
        oldzoomCommand=varargin{1};
        % Otherwise, the argument is assumed to be a oldzoom factor.
    else
        scale_factor=varargin{1};
        oldzoomCommand='scale';
    end 
    fig=get(0,'currentfigure');
    if isempty(fig)
        return
    end
    % two arg in
case 2
    fig=varargin{1};
    if ~ishandle(fig)
        error('First argument must be a figure handle.')
    end
    oldzoomCommand=varargin{2};
    % too many args
otherwise
    error(nargchk(0, 2, nargin));
end % switch nargin

%------------------------------------------------------------------------------%

oldzoomCommand=lower(oldzoomCommand);

%
% oldzoomCommand 'off'
%
if strcmp(oldzoomCommand,'off')
    % turn oldzoom off
    doZoomOff(fig);
    scribefiglisten(fig,'off');
    state = getappdata(fig,'OLDZOOMFigureState');
    if ~isempty(state)
        % since we didn't set the pointer,
        % make sure it does not get reset
        ptr = get(fig,'pointer');
        % restore figure and non-uicontrol children
        % don't restore uicontrols because they were restored
        % already when oldzoom was turned on
        uirestore(state,'nouicontrols');
        set(fig,'pointer',ptr)
        if isappdata(fig,'OLDZOOMFigureState')
            rmappdata(fig,'OLDZOOMFigureState');
        end
        % get rid of on state appdata if it exists
        % the non-existance of this appdata
        % indicates that oldzoom is off.
        if isappdata(fig,'ZoomOnState')
            rmappdata(fig,'ZoomOnState');
        end
        
    end
    % done, go home.
    return
end

%---------------------------------------------------------------------------------%

%
% set some things we need for other oldzoomCommands
%
ax=get(fig,'currentaxes');
rbbox_mode = 0;
% initialize unconstrained state
oldzoomx = 1; oldzoomy = 1;
% catch 3d oldzoom
if ~isempty(ax) & any(get(ax,'view')~=[0 90]) ...
        & ~(strcmp(oldzoomCommand,'scale')| ...
        strcmp(oldzoomCommand,'fill'))
    fZoom3d = 1;
else
    fZoom3d = 0;
end

%----------------------------------------------------------------------------------%

%
% the oldzoomCommand is 'toggle'
%
if strcmp(oldzoomCommand,'toggle'),
    state = getappdata(fig,'OLDZOOMFigureState');
    if isempty(state)
        oldzoom(fig,'on');
    else
        oldzoom(fig,'off');
    end
    % done, go home
    return
end % if

%----------------------------------------------------------------------------------%

%
% Set/check a few more things before switching to one of remaining oldzoom actions
%

% Set oldzoomx,oldzoomy and oldzoomCommand for constrained oldzooms
if strcmp(oldzoomCommand,'xdown'),
    oldzoomy = 0; oldzoomCommand = 'down'; % Constrain y
elseif strcmp(oldzoomCommand,'ydown')
    oldzoomx = 0; oldzoomCommand = 'down'; % Constrain x
end

% Catch bad argin/argout match
if (nargout ~= 0) & ...
        ~isequal(oldzoomCommand,'getmode') & ...
        ~isequal(oldzoomCommand,'getlimits') & ...
        ~isequal(oldzoomCommand,'getconnect')
    error(['OLDZOOM only returns an output if the command is getmode,' ...
            ' getlimits, or getconnect']);
end

%----------------------------------------------------------------------------------%

%
% Switch for rest of oldzoomCommands
%

switch oldzoomCommand
case 'down'
    % Activate axis that is clicked in
    allAxes = findall(datachildren(fig),'flat','type','axes');
    OLDZOOM_found = 0;
    
    % this test may be causing failures for 3d axes
    for i=1:length(allAxes)
        ax=allAxes(i);
        OLDZOOM_Pt1 = get(ax,'CurrentPoint');
        xlim = get(ax,'xlim');
        ylim = get(ax,'ylim');
        if (xlim(1) <= OLDZOOM_Pt1(1,1) & OLDZOOM_Pt1(1,1) <= xlim(2) & ...
                ylim(1) <= OLDZOOM_Pt1(1,2) & OLDZOOM_Pt1(1,2) <= ylim(2))
            OLDZOOM_found = 1;
            set(fig,'currentaxes',ax);
            break
        end
    end
    
    if OLDZOOM_found==0
        return
    end
    
    % Check for selection type
    selection_type = get(fig,'SelectionType');
    oldzoomMode = getappdata(fig,'OLDZOOMFigureMode');
    
    axz = get(ax,'ZLabel');
    
    if fZoom3d
        viewData = getappdata(axz,'OLDZOOMAxesView');
        if isempty(viewData)
            viewProps = { 'CameraTarget'...
                    'CameraTargetMode'...
                    'CameraViewAngle'...
                    'CameraViewAngleMode'};
            setappdata(axz,'OLDZOOMAxesViewProps', viewProps);
            setappdata(axz,'OLDZOOMAxesView', get(ax,viewProps));
        end
        if isempty(oldzoomMode) | strcmp(oldzoomMode,'in');
            oldzoomLeftFactor = 1.5;
            oldzoomRightFactor = .75;         
        elseif strcmp(oldzoomMode,'out');
            oldzoomLeftFactor = .75;
            oldzoomRightFactor = 1.5;
        end
        switch selection_type
        case 'open'
            set(ax,getappdata(axz,'OLDZOOMAxesViewProps'),...
                getappdata(axz,'OLDZOOMAxesView'));
        case 'normal'
            newTarget = mean(get(ax,'CurrentPoint'),1);
            set(ax,'CameraTarget',newTarget);
            camoldzoom(ax,oldzoomLeftFactor);
        otherwise
            newTarget = mean(get(ax,'CurrentPoint'),1);
            set(ax,'CameraTarget',newTarget);
            camoldzoom(ax,oldzoomRightFactor);
        end
        return
    end
    
    if isempty(oldzoomMode) | strcmp(oldzoomMode,'in');
        switch selection_type
        case 'normal'
            % Zoom in
            m = 1;
            scale_factor = 2; % the default oldzooming factor
        case 'open'
            % Zoom all the way out
            oldzoom(fig,'out');
            return;
        otherwise
            % Zoom partially out
            m = -1;
            scale_factor = 2;
        end
    elseif strcmp(oldzoomMode,'out')
        switch selection_type
        case 'normal'
            % Zoom partially out
            m = -1;
            scale_factor = 2;
        case 'open'
            % Zoom all the way out
            oldzoom(fig,'out');
            return;
        otherwise
            % Zoom in
            m = 1;
            scale_factor = 2; % the default oldzooming factor
        end
    else % unrecognized oldzoomMode
        return
    end
    
    OLDZOOM_Pt1 = get_currentpoint(ax);
    OLDZOOM_Pt2 = OLDZOOM_Pt1;
    center = OLDZOOM_Pt1;
    
    if (m == 1)
        % Zoom in
        units = get(fig,'units'); set(fig,'units','pixels')
        
        rbbox([get(fig,'currentpoint') 0 0],get(fig,'currentpoint'),fig);
        
        OLDZOOM_Pt2 = get_currentpoint(ax);
        set(fig,'units',units)
        
        % Note the currentpoint is set by having a non-trivial up function.
        if min(abs(OLDZOOM_Pt1-OLDZOOM_Pt2)) >= ...
                min(.01*[diff(get_xlim(ax)) diff(get_ylim(ax))]),
            % determine axis from rbbox 
            a = [OLDZOOM_Pt1;OLDZOOM_Pt2]; a = [min(a);max(a)];
            
            % Undo the effect of get_currentpoint for log axes
            if strcmp(get(ax,'XScale'),'log'),
                a(1:2) = 10.^a(1:2);
            end
            if strcmp(get(ax,'YScale'),'log'),
                a(3:4) = 10.^a(3:4);
            end
            rbbox_mode = 1;
        end
    end
    limits = oldzoom(fig,'getlimits');
    
case 'scale',
    if all(get(ax,'view')==[0 90]), % 2D oldzooming with scale_factor
        
        % Activate axis that is clicked in
        OLDZOOM_found = 0;
        ax = gca;
        xlim = get(ax,'xlim');
        ylim = get(ax,'ylim');
        OLDZOOM_Pt1 = [sum(xlim)/2 sum(ylim)/2];
        OLDZOOM_Pt2 = OLDZOOM_Pt1;
        center = OLDZOOM_Pt1;
        
        if (xlim(1) <= OLDZOOM_Pt1(1,1) & OLDZOOM_Pt1(1,1) <= xlim(2) & ...
                ylim(1) <= OLDZOOM_Pt1(1,2) & OLDZOOM_Pt1(1,2) <= ylim(2))
            OLDZOOM_found = 1;
        end % if
        
        if OLDZOOM_found==0, return, end
        
        if (scale_factor >= 1)
            m = 1;
        else
            m = -1;
        end
        
    else % 3D
        old_CameraViewAngle = get(ax,'CameraViewAngle')*pi/360;
        ncva = atan(tan(old_CameraViewAngle)*(1/scale_factor))*360/pi;
        set(ax,'CameraViewAngle',ncva);
        return;
    end
    
    limits = oldzoom(fig,'getlimits');
    
case 'getmode'
    state = getappdata(fig,'OLDZOOMFigureState');
    if isempty(state)
        out = 'off';
    else
        mode = getappdata(fig,'OLDZOOMFigureMode');
        if isempty(mode)
            out = 'on';
        else
            out = mode;
        end
    end
    return
    
    
case 'on',
    state = getappdata(fig,'OLDZOOMFigureState');
    if isempty(state),
        % turn off all other interactive modes
        state = uiclearmode(fig,'docontext','oldzoom',fig,'off');
        % restore button down functions for uicontrol children of the figure
        uirestore(state,'uicontrols');
        setappdata(fig,'OLDZOOMFigureState',state);
    end
    
    set(fig,'windowbuttondownfcn','oldzoom(gcbf,''down'')', ...
        'windowbuttonupfcn','ones;', ...
        'windowbuttonmotionfcn','', ...
        'buttondownfcn','', ...
        'interruptible','on');
    set(ax,'interruptible','on');
    % set an appdata so it will always be possible to 
    % determine whether or not oldzoom is on and in what
    % type of 'on' state it is.
    % this appdata will not exist when oldzoom is off
    setappdata(fig,'ZoomOnState','on');
    
    scribefiglisten(fig,'on');
    
    doZoomIn(fig)
    return
    
case 'inmode'
    oldzoom(fig,'on');
    doZoomIn(fig)
    return   
    
case 'outmode'
    oldzoom(fig,'on');
    doZoomOut(fig)
    return
    
case 'reset',
    axz = get(ax,'ZLabel');
    if isappdata(axz,'OLDZOOMAxesData')
        rmappdata(axz,'OLDZOOMAxesData');
    end
    return
    
case 'xon',
    oldzoom(fig,'on') % Set up userprop
    set(fig,'windowbuttondownfcn','oldzoom(gcbf,''xdown'')', ...
        'windowbuttonupfcn','ones;', ...
        'windowbuttonmotionfcn','','buttondownfcn','',...
        'interruptible','on');
    set(ax,'interruptible','on')
    % set an appdata so it will always be possible to 
    % determine whether or not oldzoom is on and in what
    % type of 'on' stat it is.
    % this appdata will not exist when oldzoom is off
    setappdata(fig,'ZoomOnState','xon');
    return
    
case 'yon',
    oldzoom(fig,'on') % Set up userprop
    set(fig,'windowbuttondownfcn','oldzoom(gcbf,''ydown'')', ...
        'windowbuttonupfcn','ones;', ...
        'windowbuttonmotionfcn','','buttondownfcn','',...
        'interruptible','on');
    set(ax,'interruptible','on')
    % set an appdata so it will always be possible to 
    % determine whether or not oldzoom is on and in what
    % type of 'on' state it is.
    % this appdata will not exist when oldzoom is off
    setappdata(fig,'ZoomOnState','yon');
    return
    
case 'out',
    limits = oldzoom(fig,'getlimits');
    center = [sum(get_xlim(ax))/2 sum(get_ylim(ax))/2];
    m = -inf; % Zoom totally out
    
case 'getlimits', % Get axis limits
    axz = get(ax,'ZLabel');
    limits = getappdata(axz,'OLDZOOMAxesData');
    % Do simple checking of userdata
    if size(limits,2)==4 & size(limits,1)<=2, 
        if all(limits(1,[1 3])<limits(1,[2 4])), 
            getlimits = 0; out = limits(1,:);
            return   % Quick return
        else
            getlimits = -1; % Don't munge data
        end
    else
        if isempty(limits)
            getlimits = 1;
        else 
            getlimits = -1;
        end
    end
    
    % If I've made it to here, we need to compute appropriate axis
    % limits.
    
    if isempty(getappdata(axz,'OLDZOOMAxesData')),
        % Use quick method if possible
        xlim = get_xlim(ax); xmin = xlim(1); xmax = xlim(2); 
        ylim = get_ylim(ax); ymin = ylim(1); ymax = ylim(2); 
        
    elseif strcmp(get(ax,'xLimMode'),'auto') & ...
            strcmp(get(ax,'yLimMode'),'auto'),
        % Use automatic limits if possible
        xlim = get_xlim(ax); xmin = xlim(1); xmax = xlim(2); 
        ylim = get_ylim(ax); ymin = ylim(1); ymax = ylim(2); 
        
    else
        % Use slow method only if someone else is using the userdata
        h = get(ax,'Children');
        xmin = inf; xmax = -inf; ymin = inf; ymax = -inf;
        for i=1:length(h),
            t = get(h(i),'Type');
            if ~strcmp(t,'text'),
                if strcmp(t,'image'), % Determine axis limits for image
                    x = get(h(i),'Xdata'); y = get(h(i),'Ydata');
                    x = [min(min(x)) max(max(x))];
                    y = [min(min(y)) max(max(y))];
                    [ma,na] = size(get(h(i),'Cdata'));
                    if na>1 
                        dx = diff(x)/(na-1);
                    else 
                        dx = 1;
                    end
                    if ma>1
                        dy = diff(y)/(ma-1);
                    else
                        dy = 1;
                    end
                    x = x + [-dx dx]/2; y = y + [-dy dy]/2;
                end
                xmin = min(xmin,min(min(x)));
                xmax = max(xmax,max(max(x)));
                ymin = min(ymin,min(min(y)));
                ymax = max(ymax,max(max(y)));
            end
        end
        
        % Use automatic limits if in use (override previous calculation)
        if strcmp(get(ax,'xLimMode'),'auto'),
            xlim = get_xlim(ax); xmin = xlim(1); xmax = xlim(2); 
        end
        if strcmp(get(ax,'yLimMode'),'auto'),
            ylim = get_ylim(ax); ymin = ylim(1); ymax = ylim(2); 
        end
    end
    
    limits = [xmin xmax ymin ymax];
    if getlimits~=-1, % Don't munge existing data.
        % Store limits OLDZOOMAxesData
        % store it with the ZLabel, so that it's cleared if the 
        % user plots again into this axis.  If that happens, this
        % state is cleared
        axz = get(ax,'ZLabel');
        setappdata(axz,'OLDZOOMAxesData',limits);
    end
    
    out = limits;
    return
    
case 'getconnect', % Get connected axes
    axz = get(ax,'ZLabel');
    limits = getappdata(axz,'OLDZOOMAxesData');
    if all(size(limits)==[2 4]), % Do simple checking
        out = limits(2,[1 2]);
    else
        out = [ax ax];
    end
    return
    
case 'fill',
    old_view = get(ax,'view');
    view(45,45);
    set(ax,'CameraViewAngleMode','auto');
    set(ax,'CameraViewAngle',get(ax,'CameraViewAngle'));
    view(old_view);
    return
    
otherwise
    error(['Unknown option: ',oldzoomCommand,'.']);
end

%---------------------------------------------------------------------------------%

%
% Actual oldzoom operation
%

if ~rbbox_mode,
    xmin = limits(1); xmax = limits(2); 
    ymin = limits(3); ymax = limits(4);
    
    if m==(-inf),
        dx = xmax-xmin;
        dy = ymax-ymin;
    else
        dx = diff(get_xlim(ax))*(scale_factor.^(-m-1)); dx = min(dx,xmax-xmin);
        dy = diff(get_ylim(ax))*(scale_factor.^(-m-1)); dy = min(dy,ymax-ymin);
    end
    
    % Limit oldzoom.
    center = max(center,[xmin ymin] + [dx dy]);
    center = min(center,[xmax ymax] - [dx dy]);
    a = [max(xmin,center(1)-dx) min(xmax,center(1)+dx) ...
            max(ymin,center(2)-dy) min(ymax,center(2)+dy)];
    
    % Check for log axes and return to linear values.
    if strcmp(get(ax,'XScale'),'log'),
        a(1:2) = 10.^a(1:2);
    end
    if strcmp(get(ax,'YScale'),'log'),
        a(3:4) = 10.^a(3:4);
    end
    
end

% Check for axis equal and update a as necessary
if strcmp(get(ax,'plotboxaspectratiomode'),'manual') & ...
        strcmp(get(ax,'dataaspectratiomode'),'manual')
    ratio = get(ax,'plotboxaspectratio')./get(ax,'dataaspectratio');
    dx = a(2)-a(1);
    dy = a(4)-a(3);
    [kmax,k] = max([dx dy]./ratio(1:2));
    if k==1
        dy = kmax*ratio(2);
        a(3:4) = mean(a(3:4))+[-dy dy]/2;
    else
        dx = kmax*ratio(1);
        a(1:2) = mean(a(1:2))+[-dx dx]/2;
    end
end

% Update circular list of connected axes
list = oldzoom(fig,'getconnect'); % Circular list of connected axes.
if oldzoomx,
    if a(1)==a(2), return, end % Short circuit if oldzoom is moot.
    set(ax,'xlim',a(1:2))
    h = list(1);
    while h ~= ax,
        set(h,'xlim',a(1:2))
        % Get next axes in the list
        hz = get(h,'ZLabel');
        next = getappdata(hz,'OLDZOOMAxesData');
        if all(size(next)==[2 4]), h = next(2,1); else h = ax; end
    end
end
if oldzoomy,
    if a(3)==a(4), return, end % Short circuit if oldzoom is moot.
    set(ax,'ylim',a(3:4))
    h = list(2);
    while h ~= ax,
        set(h,'ylim',a(3:4))
        % Get next axes in the list
        hz = get(h,'ZLabel');
        next = getappdata(hz,'OLDZOOMAxesData');
        if all(size(next)==[2 4]), h = next(2,2); else h = ax; end
    end
end

%----------------------------------------------------------------------------------%

%
% Helper Functions
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function p = get_currentpoint(ax)
%GET_CURRENTPOINT Return equivalent linear scale current point
p = get(ax,'currentpoint'); p = p(1,1:2);
if strcmp(get(ax,'XScale'),'log'),
    p(1) = log10(p(1));
end
if strcmp(get(ax,'YScale'),'log'),
    p(2) = log10(p(2));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function xlim = get_xlim(ax)
%GET_XLIM Return equivalent linear scale xlim
xlim = get(ax,'xlim');
if strcmp(get(ax,'XScale'),'log'),
    xlim = log10(xlim);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ylim = get_ylim(ax)
%GET_YLIM Return equivalent linear scale ylim
ylim = get(ax,'ylim');
if strcmp(get(ax,'YScale'),'log'),
    ylim = log10(ylim);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function doZoomIn(fig)
set(findall(fig,'Tag','figToolZoomIn'),'State','on');   
set(findall(fig,'Tag','figToolZoomOut'),'State','off');      
setappdata(fig,'OLDZOOMFigureMode','in');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function doZoomOut(fig)
set(findall(fig,'Tag','figToolZoomIn'),'State','off');   
set(findall(fig,'Tag','figToolZoomOut'),'State','on');      
setappdata(fig,'OLDZOOMFigureMode','out');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function doZoomOff(fig)
set(findall(fig,'Tag','figToolZoomIn'),'State','off');
set(findall(fig,'Tag','figToolZoomOut'),'State','off');   
if ~isempty(getappdata(fig,'ZoomFigureMode'))
    rmappdata(fig,'OLDZOOMFigureMode');
end
