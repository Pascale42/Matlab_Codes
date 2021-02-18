function out = VonMisesFit(varargin)
% out = VonMisesFitCorr(th) or VonMisesFitCorr(r,n)
% fits a Von Mises distribution by conditional  maximum likelihood
% this solves the problem of bias of r for different sample size
% see Schou , Ekonometrika 1978


if nargin==1
    th = varargin{1};
    [mu r] = circmean(th);
    n = length(th);
else
    r= varargin{1};
    n = varargin{2};
end
%nnow the usual ML estimate
kml = BesselRatInv(r);

if n<=15
    if kml<2
        out.kml = max(kml -2/(n*kml),0);
    else
        out.kml = (n-1)^3*kml/(n^3+n);
    end
 else
    out.kml = kml;
end

out.knr = NewtonRaphson(out.kml/5,r,n);
out.kfz = fzero(@(x) SchouFun(x,r,n),out.kml/5);
%keyboard
U = inline('R.*vmA(k*R)./vmA(k)-n','k','R','n'); %g(k,R,n)
H = inline('R./k.*vmA(k*R)./vmA(k).*(vmB(k)-vmB(k*R))','k','R'); % dg/dk

f = @(x) U(x,r,n);
df = @(x) H(x,r);
[out.knr1,y0,err]=newton(f,df,out.kml);

out.r=r;

return


%function x =NewtonRaphson(xo,r,n M, MaxTol)
function x =NewtonRaphson(xo,r,n,varargin)
[M, MaxTol] = DefaultArgs(varargin,{50,1e-8});
R=r*n;

%for the ML estimate
%U = inline('r*n-n*bessel(1,k)/bessel(0,k)','k','r','n');
%H = inline('n*(bessel(1,k)/bessel(0,k).^2+bessel(1,k)/bessel(0,k)/k-1)','k','n');


U = inline('R.*vmA(k*R)./vmA(k)-n','k','R','n'); %g(k,R,n)
H = inline('R./k.*vmA(k*R)./vmA(k).*(vmB(k)-vmB(k*R))','k','R'); % dg/dk

i=0;
xcur=xo;
while i<M & abs(U(xcur,R,n))>MaxTol
   % fprintf('iteration # %d k = %2.5f\n',i,xcur);
    i=i+1;
    xprev = xcur;
    xcur = xprev-U(xprev,R,n)./H(xprev,R);
   
end
if abs(U(xcur,R,n))<=MaxTol
    x = xcur;
else
    warning('algorythm exceeded max number of iterations');
    x=xcur;
    return;
end
%return

function a = A(k)
a  = bessel(1,k)/bessel(0,k);
return
