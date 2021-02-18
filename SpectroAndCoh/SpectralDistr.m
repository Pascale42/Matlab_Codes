function [yh, by, phih, bphi] = SpectralDistr(y,f,t,varargin)
%function [yh, by, phih, bphi] = SpectralDistr(y,f,t,phi,norm,period, nbins)
% period - in seconds
[phi,norm,period,nbins] = DefaultArgs(varargin,{[],0,[],100});
nCh = size(y,4);
nF = length(f);
nT=length(t);
yh =[]; phih = [];
for Ch1=1:nCh
    for Ch2=Ch1:nCh
        p=squeeze(y(:,:,Ch1,Ch2));
        if ~isempty(period)
            mask = repmat(WithinRanges(t,period),1,nF);
            nsz = [sum(mask(:,1)) size(p,2)];
            p = p(logical(mask));
            p = reshape(p,nsz);
        end
        
%        p = OutlierRemoveAuto(p,0.1);
        thr = prctile(p(:),[1 99.9]);
        p(p(:)>thr(2)) = thr(2);
        p(p(:)<thr(1)) = thr(1);
        if Ch1==Ch2
            p = log(p);
        else
            p = atanh(p);
        end
        [h,b] = hist(p,nbins);
        if norm
        h=h./repmat(sum(h),nbins,1);
        end
        if nargout<1
            figure(2212)
            if norm
            nh=h./(repmat(sum(h,2),1,nF)+eps);
            end
            subplot(nCh,nCh,(Ch1-1)*nCh+Ch2);
            if norm
                imagesc(f,b,nh);axis xy
            else
                imagesc(f,b,h);axis xy
            end
%             if Ch1~=Ch2
%                 tit  = ['Coherence ' num2str(Ch1) ' - ' num2str(Ch2)];
%                 ax = str2num(get(gca,'YTickLabel'));
%                 set(gca,'YTickLabel',num2str(tanh(ax)));
%             else
%                 tit  = ['Power ' num2str(Ch2)];
% 
%             end
            %title(tit);
        else
            yh(:,:,Ch1,Ch2) = h;
            by(:,Ch1,Ch2) = b;
        end

        %now for the phase
        if ~isempty(phi)
            p=squeeze(phi(:,:,Ch1,Ch2));
            if ~isempty(period)
                mask = repmat(WithinRanges(t,period),1,nF);
                nsz = [sum(mask(:,1)) size(p,2)];
                p = p(logical(mask));
                p = reshape(p,nsz);

            end
            if Ch1==Ch2
                continue;
            end
            [h,b] = hist([p; p+2*pi],nbins*2);
            h=h./repmat(sum(h),nbins*2,1);
            if nargout<1
               % figure(2212)
                nh=h./(repmat(sum(h,2),1,nF)+eps);
                subplot(nCh,nCh,(Ch2-1)*nCh+Ch1);
                if norm
                    imagesc(f,b,nh);axis xy
                else
                    imagesc(f,b,h);axis xy
                end
                tit  = ['phase ' num2str(Ch1) ' - ' num2str(Ch2)];
                title(tit);
            else
                phih(:,:,Ch1,Ch2) = h;
                bphi(:,Ch1,Ch2) = b*180/pi;
            end


        end
    end
end