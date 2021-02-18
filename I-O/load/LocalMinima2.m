%function Mins = LocalMinima2(p,AmpThr,NotCloserThen)
% gives the coordinates of local minima of surface given by matrix p
% 
function Mins = LocalMinima2(p,AmpThr,NotCloserThan)

if length(NotCloserThan)==1 NotCloserThan = [1 1]*NotCloserThan;end

%here we will try to do left, right,up and down shifts to get the local
%minima : right and down are positive , left and up are negative
[nx,ny] = size(p);
uCond = [ones(1,ny); diff(p,1,1)]<0;
dCond = [diff(p,1,1); -ones(1,ny)]>0;

lCond = [ones(nx,1) diff(p,1,2)]<0;
rCond = [diff(p,1,2) -ones(nx,1)]>0;

ampCond = p<AmpThr;

[Mins(:,1), Mins(:,2)] = find(dCond & uCond & rCond & lCond & ampCond);

while 1
    %cheap approximate way
    nMins= size(Mins,1);
    Disti = abs(repmat(Mins(:,1),1,nMins)-repmat(Mins(:,1)',nMins,1));
    Distj = abs(repmat(Mins(:,2),1,nMins)-repmat(Mins(:,2)',nMins,1));
    [TooClose1 TooClose2] = find(tril(Disti<NotCloserThan(1) & Distj<NotCloserThan(2),-1));
    nTooClose=length(TooClose1);
    if isempty(TooClose1)
        break;
    end
    ii1=sub2ind([nx ny],Mins(TooClose1,1),Mins(TooClose1,2));
    ii2=sub2ind([nx ny],Mins(TooClose2,1),Mins(TooClose2,2));
    Vals = [p(ii1) p(ii2)];
    [dummy Offset] = max(Vals,[],2);
    TooClose = [TooClose1(:) TooClose2(:)];
    Delete = TooClose(sub2ind([nTooClose,2],[1:nTooClose]',Offset));
    Mins(unique(Delete),:) = [];
end

if nargout<1
    imagesc(p');axis xy
    hold on
    plot(Mins(:,1),Mins(:,2),'ow','MarkerSize',15)
    plot(Mins(:,1),Mins(:,2),'kx','MarkerSize',15)
end



return

% compute derivatives in x- and y-directions
[dp_dx,dp_dy] = gradient(p,x,y);

% compute second derivatives
[d2p_dx2,d2p_dy_dx] = gradient(dp_dx,x,y);
[d2p_dx_dy,d2p_dy2] = gradient(dp_dy,x,y);

%lines of zero derivative 
cx = contourc(x,y,dp_dx,[0 0]);
cy = contourc(x,y,dp_dy,[0 0]);

%now go through curves of zero derivatives and find their intersections

%get actual curves from crazy contourc output
curvex={}; curvey={};
cnt=1; headptr =1;
while headptr<size(cx,2)
    clen = cx(2,headptr);
    curvex{cnt} = cx(:,headptr+[1:clen]);
    headptr = headptr+clen+1;
    cnt=cnt+1;
end

cnt=1; headptr =1;
while headptr<size(cy,2)
    clen = cy(2,headptr);
    curvey{cnt} = cy(:,headptr+[1:clen]);
    headptr = headptr+clen+1;
    cnt=cnt+1;
end

%now find intersections

intersects=[];
for i=1:length(curvex)
    for j=1:length(curvey)
 %       clf
%        plot(curvex{i}(1,:),curvex{i}(2,:),'b');hold on;         
  %      plot(curvey{j}(1,:),curvey{j}(2,:),'g');
      %  [ix iy] = curveintersect(curvex{i}(1,:),curvex{i}(2,:),curvey{j}(1,:),curvey{j}(2,:));
        intcoord = CurvesIntersection(curvex{i}',curvey{j}');
   %     plot(intcoord(:,1),intcoord(:,2),'rx','MarkerSize',10);
    %    pause
        intersects=[intersects; intcoord];

    end
end
extrem=[];
%[dummy extrem(:,1)] = NearestNeighbour(x,intersects(:,1));
%[dummy extrem(:,2)] = NearestNeighbour(y,intersects(:,2));

%interpolate points amplitudes and derivatives
amp = interp2(x,y,p,intersects(:,1),intersects(:,2));
d2px = interp2(x,y,d2p_dx2,intersects(:,1),intersects(:,2));
d2py = interp2(x,y,d2p_dy2,intersects(:,1),intersects(:,2));

mins = intersects(d2px>0 & d2py>0 & amp<AmpThr(1),:);

maxs = intersects(d2px<0 & d2py<0 & amp>AmpThr(2),:);

%different signs - saddle , but we don't care for that now

if nargout<1
    imagesc(x,y,p);axis xy
    hold on
    plot(mins(:,1),mins(:,2),'xr','MarkerSize',15)
    plot(maxs(:,1),maxs(:,2),'xb','MarkerSize',15)
end


return

% Let's plot those regions using yellow for maxima, gray for saddle points
% and cyan for minima and then overlay the contour lines from the first
% figure

% % create array with maxima and minima
p = peaks(101);
x = 0:100; % x values
y = 0:100; % y values

LocalMinima2(p,x,y,[0 0]);