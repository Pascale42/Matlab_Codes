%function [pvm, pu, r] = VonMisesRPdf(k,n,pu)
function [pvm, pu, r] = VonMisesRPdf(k,n,varargin)
[pu] = DefaultArgs(varargin,{[]});

r=[0:0.001:1];
R=r.*n;
%y= besseli(0,k.*R)./(besseli(0,k).^n).*R;
y= besseli(0,k.*R)./(besseli(0,k).^n);
% y= y*UniformRPdf
% we need to express analityically, or estimate with quadratures or
% montecarlo the UniformRPdf

%yr =hist(Ru,R);
%yr = yr./sum(yr);
if isempty(pu)
    thu = rand(n,100000)*2*pi;%VonMisesRnd(zeros(n,1000),k*ones(n,1000));
    Ru = abs(sum(exp(i*thu)));
    [pu,xi]=ksdensity(Ru,R,'support','positive');
    pu=pu./sum(pu);
    pu =pu(:);
elseif length(pu)==1
    pu=[];
    for Rcur=R
        pu(end+1,1)= Rcur.*quad(hnfun, 0, 1e8);
      %  pR(end+1,1) = bessel(0,k*R)*hn./(bessel(0,k).^n);
    end
    
end
pvm = y(:).*pu;
%pu=pu(:); pvm=pvm(:);
if nargout<1
    figure
    plot(r,[pu pvm]);
    legend('uniform','von Mises');
    title(sprintf('k=%f, n=%d',k,n));    
end


    function y= hnfun(u)
        y = @(u) u.*bessel(0,u.*Rcur).*(bessel(0,u).^n);
    end

end