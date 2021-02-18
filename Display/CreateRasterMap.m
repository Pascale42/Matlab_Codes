%function Map = CreateRasterMap(action,GridDims,Name)
% creates rater map (for the cursosr or whatever)
% default action is start (can leave empty), 
% GridDims - x and y dimmensions of the grid (16 by 16 for cursor)
% left click - 1 (black), right click -2 (white),
% the rest - NaN. bressing b - you can right click and put news NaNs
% q - quit and save the results in the Name.map ASCII file
function Map = CreateRasterMap(varargin)
global gCreateRasterMap;
[action,GridDims,Name] = DefaultArgs(varargin,{'start',[16 16],'Default'});

switch action 
case 'start'
	gCreateRasterMap = struct;
	gCreateRasterMap.figh = figure;
	axis([0 GridDims(1) 0 GridDims(2)]);
	%init matrix
    if FileExists([Name '.map']);
        map=load([Name '.map']);
        gCreateRasterMap.Map = flipud(map)';
    else
        gCreateRasterMap.Map = NaN*ones(GridDims);
    end
	gCreateRasterMap.DoBorders = 0;
	gCreateRasterMap.Name =Name;
	CreateRasterMap('update');
	
	% if click 
	set(gcf,'WindowButtonDownFcn','CreateRasterMap(''click'')');
	% if release
	set(gcf,'WindowButtonUpFcn','CreateRasterMap(''unclick'')');
	% if press key
	set(gcf,'KeyPressFcn','CreateRasterMap(''key'')');

case 'update'
	imagesc([1:GridDims(1)],[1:GridDims(2)],gCreateRasterMap.Map');
	caxis([0 2]);
	axis xy
	set(gca,'XTick',[1-0.5:1:GridDims(1)+0.5]);
	set(gca,'YTick',[1-0.5:1:GridDims(1)+0.5]);
	grid on
	set(gca,'GridLineStyle','-');

case 'move'
	whatbutton = get(gcf,'SelectionType');
	mousecoord = get(gca,'CurrentPoint');
	coord = round(mousecoord(1,1:2));
	coord = min(coord,GridDims);
	coord = max(coord,[1 1]);
	%coord = [floor(realcoord)-
	switch whatbutton
	case 'normal' %left  add point
		if ~gCreateRasterMap.DoBorders
			gCreateRasterMap.Map(coord(1),coord(2)) = 1;		
		else
			gCreateRasterMap.Map(coord(1),coord(2))=gCreateRasterMap.OldMap(coord(1),coord(2));
		end
	case 'alt' % right - erase point
		if ~gCreateRasterMap.DoBorders
			gCreateRasterMap.Map(coord(1),coord(2)) = 2;		
		else
			gCreateRasterMap.Map(coord(1),coord(2)) = NaN;
		end		
	end	
	 CreateRasterMap('update');
case 'click'
	whatbutton = get(gcf,'SelectionType');
	mousecoord = get(gca,'CurrentPoint');
	coord = round(mousecoord(1,1:2));
	coord = min(coord,GridDims);
	coord = max(coord,[1 1]);
	switch whatbutton
	case 'normal' %left  add point
		if ~gCreateRasterMap.DoBorders
			gCreateRasterMap.Map(coord(1),coord(2)) = 1;		
		else
			gCreateRasterMap.Map(coord(1),coord(2))=gCreateRasterMap.OldMap(coord(1),coord(2));
		end
	case 'alt' % right - erase point
		if ~gCreateRasterMap.DoBorders
			gCreateRasterMap.Map(coord(1),coord(2)) = 2;		
		else
			gCreateRasterMap.Map(coord(1),coord(2)) = NaN;
		end		
	%case 'open' % double click - save and quit
	%case 'extend' %middle click - erase all
	end
	CreateRasterMap('update');
	set(gcf,'WindowButtonMotionFcn','CreateRasterMap(''move'')');
	
case 'unclick'
	set(gcf,'WindowButtonMotionFcn','');
	
case 'key'
	whatkey = get(gcf,'CurrentCharacter');
	switch whatkey
	case 'q'
		%global Map
		msave([gCreateRasterMap.Name '.map'],flipud(gCreateRasterMap.Map'));
		close(gCreateRasterMap.figh);
		figure
		imagesc(gCreateRasterMap.Map');caxis([0 2]);
		axis xy
		title('Here is your raster map');
		set(gcf,'Pointer','custom');
		set(gcf,'PointerShapeCData',flipud(gCreateRasterMap.Map'));
		clear gCreateRasterMap;
		%varargout{1} = gCreateRasterMap.Map;
		return
	case 'e' % erase all
		gCreateRasterMap.Map = NaN*ones(GridDims);
		CreateRasterMap('update');
	case 'b' % make borders
		gCreateRasterMap.DoBorders = 1;
		gCreateRasterMap.OldMap = gCreateRasterMap.Map;
	end

%case 'return'
	
end
