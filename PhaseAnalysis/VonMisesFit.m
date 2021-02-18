% [mu, k, muerr] = VonMisesFit(th)
% fits a Von Mises distribution by maximum likelihood
%
% muerr WILL BE STANDARD ERROR AS SOON AS I DO IT

function [mu, k, muerr, kerr] = VonMisesFit(th)

[mu r] = circmean(th);

kml = BesselRatInv(r);

n = length(th(:));
if n<=15
    if kml<2
        k = max(kml -2/(n*kml),0);
    else
        k = (n-1)^3*kml/(n^3+n);
    end
else
   k = kml;
  % k = NewtonRaphson(kml,r,n);
  % fprintf('k0 = %2.2f, k = %2.2f\n', k0,k);
end

% muerr = 1./n./k/bessel(1,k)*bessel(0,k);
% kerr = 1./(n*(1-bessel(1,k)/bessel(0,k).^2-bessel(1,k)/bessel(0,k)./k));
%keyboard

return

% function a = A(k)
% a  = bessel(1,k)/bessel(0,k);
% return

%function x =NewtonRaphson(xo,r,n M, MaxTol)
function x =NewtonRaphson(xo,r,n,varargin)
[M, MaxTol] = DefaultArgs(varargin,{50,1e-8});


U = inline('r*n-n*bessel(1,k)/bessel(0,k)','k','r','n');
H = inline('n*(bessel(1,k)/bessel(0,k).^2+bessel(1,k)/bessel(0,k)/k-1)','k','n');

i=0;
xcur=xo;
while i<M & abs(U(xcur,r,n))>MaxTol
    fprintf('iteration # %d k = %2.5f\n',i,xcur);
    i=i+1;
    xprev = xcur;
    xcur = xprev-U(xprev,r,n)./H(xprev,n);
   
end
if abs(U(xcur,r,n))<=MaxTol
    x = xcur;
else
    error('algorythm exceeded max number of iterations');
end
return

