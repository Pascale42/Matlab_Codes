% function PhPowCorr(per)
% % per = nt x nf x ntapers x 2
% 
% [nt,nf,ntap,nch] = size(per);
% csd = sq(mean(mper(:,:,:,1).*conj(mper(:,:,:,2)),3));
% pow(:,:,1) = sq(abs(mean(csd(:,:,:,1).^2),3));
% pow(:,:,2) = sq(abs(mean(csd(:,:,:,2).^2),3));
% 
% ph = angle(csd);
% 
% mph = angle(mean(exp(sqrt(-1)*ph)));
% mR = abs(mean(exp(sqrt(-1)*ph)));
% 
% cph = mod(dotdot(ph,'-',mph)+pi,2*pi)-pi;
% 
% [hst cphi] = histcI(1-cos(cph),linspace(-pi,pi,7));
% 
% for l=1:nf
%     for k=1:6
%        myamp = amp(cphi(:,l)==k);
%        rcorr(l,k)= RankCorrelation(, 
%     end