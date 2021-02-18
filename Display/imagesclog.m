% function zlog=imagesclog(x,y,z,plotit);
%
% x and y are vectors
% z is a matrix of size [length(x) length(y)]
% plotit = 1; plots directly
% plot output : set(gca,'XScale','log')

function zlog=imagesclog(x,y,z,plotit)

xlog=exp(linspace(log10(x(1)),log(x(end)),length(x)))';

r=length(x); c=length(y);

zlog=zeros(r,c);


for n=1:c
  zlog(:,n)=interp1(x,z(:,n),xlog,'linear');
end

if plotit ==1

figure('Name','Image log', 'NumberTitle','off');

imagesc(x,y,zlog'); set(gca,'XScale','log')

end
