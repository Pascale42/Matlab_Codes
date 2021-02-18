function [x,y]=curveintersect(varargin)
% Curve Intersections.
% [X,Y]=CURVEINTERSECT(H1,H2) or [X,Y]=CURVEINTERSECT([H1 H2]) finds the
% intersection points of the two curves on the X-Y plane identified
% by the line or lineseries object handles H1 and H2.
%
% [X,Y]=CURVEINTERSECT(X1,Y1,X2,Y2) finds the intersection points of the
% two curves described by the vector data pairs (X1,Y1) and (X2,Y2).
%
% X and Y are empty if no intersection exists.
%
% PLOT(X,Y,'o') marks the intersections with circles.

% D.C. Hanselman, University of Maine, Orono, ME 04469
% Mastering MATLAB 7
% 2005-01-06

[x1,y1,x2,y2]=local_parseinputs(varargin{:});

if ~isequal(x1,x2)
   xx=unique([x1 x2]); % get unique data points
   xx=xx(xx>=max(min(x1),min(x2)) & xx<=min(max(x1),max(x2)));
   if numel(xx)<2
      x=[];
      y=[];
      return
   end
   x1 = makeunique(x1);   x2 = makeunique(x2);
   y1 = makeunique(y1);   y2 = makeunique(y2);
   yy=interp1(x1,y1,xx)-interp1(x2,y2,xx);
   
      
else
   xx=x1;
   yy=y1-y2;
end
x=mminvinterp(xx,yy,0); % find zero crossings of difference
if ~isempty(x)
   y=interp1(x1,y1,x);
else
   x=[];
   y=[];
end

%--------------------------------------------------------------------------
function ux = makeunique(x)
%makes all elements unique by adding random jitter to repeats
ux=x(:);
[u,i,j] = unique(x);
if length(u)==length(x) return; end
repeats = setdiff(j,i);
minDiff = min(abs(diff(u)));
ux(repeats) = ux(repeats)+(rand(length(repeats),1)-0.5)*0.001*minDiff;

%--------------------------------------------------------------------------
function [xo,yo]=mminvinterp(x,y,yo)
%MMINVINTERP 1-D Inverse Interpolation. From the text "Mastering MATLAB 7"
% [Xo, Yo]=MMINVINTERP(X,Y,Yo) linearly interpolates the vector Y to find
% the scalar value Yo and returns all corresponding values Xo interpolated
% from the X vector. Xo is empty if no crossings are found. For
% convenience, the output Yo is simply the scalar input Yo replicated so
% that size(Xo)=size(Yo).
% If Y maps uniquely into X, use INTERP1(Y,X,Yo) instead.
%
% See also INTERP1.

if nargin~=3
   error('Three Input Arguments Required.')
end
n = numel(y);
if ~isequal(n,numel(x))
   error('X and Y Must have the Same Number of Elements.')
end
if ~isscalar(yo)
   error('Yo Must be a Scalar.')
end

x=x(:); % stretch input vectors into column vectors
y=y(:);

if yo<min(y) || yo>max(y) % quick exit if no values exist
   xo = [];
   yo = [];
else                      % find the desired points
   
   below = y<yo;          % True where below yo 
   above = y>=yo;         % True where at or above yo
   
   kth = (below(1:n-1)&above(2:n))|(above(1:n-1)&below(2:n)); % point k
   kp1 = [false; kth];                                        % point k+1
   
   alpha = (yo - y(kth))./(y(kp1)-y(kth));% distance between x(k+1) and x(k)
   xo = alpha.*(x(kp1)-x(kth)) + x(kth);  % linearly interpolate using alpha
   
   yo = repmat(yo,size(xo)); % duplicate yo to match xo points found
end 
%--------------------------------------------------------------------------
function [x1,y1,x2,y2]=local_parseinputs(varargin)

if nargin==1 % [X,Y]=CURVEINTERSECT([H1 H2])
   arg=varargin{1};
   if numel(arg)==2 && ...
      all(ishandle(arg)) && all(strcmp(get(arg,'type'),'line'))
      data=get(arg,{'XData','YData'});
      [x1,x2,y1,y2]=deal(data{:});
   else
      error('Input Must Contain Two Handles to Line Objects.')
   end
elseif nargin==2 % [X,Y]=CURVEINTERSECT(H1,H2)
   arg1=varargin{1};
   arg2=varargin{2};
   if numel(arg1)==1 && ishandle(arg1) && strcmp(get(arg1,'type'),'line')...
   && numel(arg2)==1 && ishandle(arg2) && strcmp(get(arg2,'type'),'line')
      
      data=get([arg1;arg2],{'XData','YData'});
      [x1,x2,y1,y2]=deal(data{:});
   else
      error('Input Must Contain Two Handles to Line Objects.')
   end
elseif nargin==4
   [x1,y1,x2,y2]=deal(varargin{:});
   if ~isequal(numel(x1),numel(y1))
      error('X1 and Y1 Must Contain the Same Number of Elements.')
   elseif ~isequal(numel(x2),numel(y2))
      error('X2 and Y2 Must Contain the Same Number of Elements.')
   end
   x1=reshape(x1,1,[]); % make data into rows
   x2=reshape(x2,1,[]);
   y1=reshape(y1,1,[]);
   y2=reshape(y2,1,[]);
else
   error('Incorrect Number of Input Arguments.')
end
if numel(x1)<2 || numel(x2)<2 || numel(y1)<2 || numel(y2)<2
   error('At Least Two Data Points are Required for Each Curve.')
end