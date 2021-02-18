%function trmom = TriganomMoments(th)
%computes various stats on the phase vector( matrix)
function out = TriganomMoments(th)
x1 = exp(i*th);
x2 = exp(2*i*th);
mx1 = mean(x1);
mx2 = mean(x2);

th1 = angle(mx1);
th2 = angle(mx2);

m1 = abs(mx1);
m2 = abs(mx2);
circvar = 1-m1;
circstd = sqrt(-2*log(m1));
circdisp = (1-m2)/2./m1.^2;
circskew = m2.*sin(th2-2*th1)./circvar.^(3/2);
circkurt = (m2.*cos(th2-2*th1)-m1.^4)./circvar.^2;

out.mean = th1;
out.var = circvar;
out.std = circstd;
out.disp = circdisp;
out.skew = circskew;
out.kurt = circkurt;


