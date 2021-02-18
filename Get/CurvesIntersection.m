%IntCoord = function CurvesIntresection(c1,c2)
function IntCoord = CurvesIntresection(c1,c2)
% searches for intersection of the curves using crazy idea of inpolygon
% function usage :))
%assumes that edges are all outside the curves closures .

%first find the boreders of the data from both curves
maxx = max([c1(:,1);c2(:,1)]);
minx = min([c1(:,1);c2(:,1)]);
maxy = max([c1(:,2);c2(:,2)]);
miny = min([c1(:,2);c2(:,2)]);
Center = [(minx+maxx)/2, (miny+maxy)/2];
Corners = [minx*0.99 minx*0.99 maxx*1.01 maxx*1.01; miny*0.99 maxy*1.01 maxy*1.01 miny*0.99]';
n1=length(c1);n2=length(c2);

if n1>n2
    Data = c1;
    Border = c2;
else
    Data = c2;
    Border = c1;
end

%now need to close the border and make it 

%now assume that curves end are ending at the borders 
%have to fix it otherwise

%find the configuration of the ends

AngleBorder = atan2(Border([1 end],2)-Center(2), Border([1 end],1)-Center(1));
%AngleData = atan2(Data([1 end],2)-Center(2), Data([1 end],1)-Center(1));
AngleCorners = atan2(Corners(:,2)-Center(2), Corners(:,1)-Center(1));

Points=[Border([1 end],:); Corners];

%order the border ends and corners
[s si] = sort(mod([AngleBorder; AngleCorners]-AngleBorder(1),2*pi));
%indexes : 1 2 -Border, 3:6 - corners (ld, lu, ru, rd)

newPoints = Points(flipud(si(1:find(si==2)-1)),:);
Border = [Border; newPoints];

In = inpolygon(Data(:,1),Data(:,2),Border(:,1),Border(:,2));
IntersInd = find(abs(diff(In)));

IntCoord = Data(IntersInd,:);
