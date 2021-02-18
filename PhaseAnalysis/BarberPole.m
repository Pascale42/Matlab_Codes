% [b k berr] = BarberPole(th, y)
%
% fits a model of the form th ~ VonMises(y*b(1) + b(2), k)
%
% berr returns a standard error estimate for b
%
% See Fisher, Statistical analysis of circular data, sec 6.4.2

function [b, k, berr] = BarberPole(th, y)

if min(size(th))>1 | min(size(y))>1
	error('inputs must be vectors')
end
th = th(:);
y = y(:);

if size(y,1)~=size(th,1)
	error('number of points do not match');
end

% test out a grid.
nGrid=30;
yr = max(y)-min(y);
ar = 0:2*pi/nGrid:2*pi;
br = (-2*pi:8*pi/nGrid:2*pi)/yr;

th0 = repmat(y*br,[1 1 length(ar)]) + repmat(permute(ar,[1 3 2]), [length(y) length(br) 1]);
dth = repmat(th, [1 length(br) length(ar)]) - th0;

l = permute(sum(cos(dth),1), [2 3 1]);

% find max val on grid (b0,a0)
[dummy ind] = max(l(:));
[bi ai] = ind2sub(size(l), ind);
b0 = br(bi);
a0 = ar(ai);


% now get serious with a search from this start point
f = inline('-sum(cos(P1-P2*x(1)-x(2)))',2);
b = fminsearch(f, [b0;a0], [], th,y);

% now find k.

dth0 = th-b(1)*y-b(2);

r = mean(cos(dth0));
k = BesselRatInv(r);

yPlus = [y, ones(length(y),1)];
E2 = yPlus'*yPlus;

berr = sqrt(diag(inv(E2))/k/r);

% plot grid
subplot(2,1,1)
hold off
imagesc(ar,br,l);

% plot fit
subplot(2,1,2);
hold off
yu = unique(y);
plot(y, mod(th+pi,2*pi)-pi, '.', yu, mod(b(2)+yu*b(1)+pi,2*pi)-pi);

drawnow

return

[nPoints nDims] = size(y);

if nargin<3
    b0 = zeros(nDims,1);
end

[m r] = circmean(th);
f = inline('-sum(cos(P1-P2*x))',2);
b = fminsearch(f, b0, [], th,y);

plot(y(:,1),th, '.', y(:,1), y*b);

return

if size(th,2)>1
	error('th must be a vector')
end

if size(x,1)~=size(th,1)
	error('number of points do not match');
end

[nPoints nDims] = size(x);

[m r] = circmean(th);
b = zeros(nDims,1);
b = [.5;0]
k = 0;
bold = inf*ones(size(b));
kold = k;

tol = 1e-7;

xt = x';
while (sum(abs(b-bold)) + abs(k-kold)) > tol
	bold = b;
	kold = k;

	g = 2*atan(x*b);

	[mu r] = circmean(th-g);
	k = BesselRatInv(r);

	bbold = inf*ones(size(b));
	while sum(abs(b-bbold))>tol;

		bbold = b;
		xb = x*b;
		g = 2*atan(xb);
        dth = th-mu-g;
        dg = 2./(1+xb.^2);
        
    %    D = x'*(sin(dth).*dg)
     %   DD = x'*sparse(1:nPoints,1:nPoints,dg.*cos(dth) + (dg.^2).*xb.*sin(dth))*x
        
      %  db = DD\D;
        
        ddg = -4*xb./(1+xb.^2).^2;
		G2 = sparse(1:nPoints,1:nPoints,dg.^2);
		G = sparse(1:nPoints,1:nPoints,dg);
        GG = sparse(1:nPoints,1:nPoints,ddg);
		y = (sin(th-mu-g)./dg)/r;

%         top = sin(th-mu-g)'*G*x;
%         bot = xt*(cos(th-mu-g)'*G + sin(th-mu-g)'*GG)*x;
 		db = (xt*G*x)\(xt*G*sin(dth)/r);

 		plot(x(:,1), th, '.', x(:,1), mu + 2*atan(xb));

		b = b+db;

%		keyboard %pause
    disp(b')
end
	drawnow
    pause
end

f = inline('-sum(cos(P1-P2*x))',2)
b1 = -2:.1:2;
b2 = -2*pi:.1:2*pi;
clear l
for i=1:length(b1)
    for j=1:length(b2)
        l(i,j) = sum(cos(th-x*[b1(i);b2(j)]));
    end
end
imagesc(b1,b2,l');
b0 = fminsearch(f, [0;0], [], th,x)