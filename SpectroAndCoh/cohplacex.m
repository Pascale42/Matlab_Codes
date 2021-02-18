fn = 'n1';
 
chans =[1:4:27];
eeg = readmulti([fn '.eeg'],31,chans);

whl = load([fn '.whl']);

nCh = length(chans);

WinLength = 64; % i.e to convert from eeg to whl , but twice bigger
nOverlap = WinLength /2; % so that we get same # of samples
[y, f, t]=mtchgband(eeg,2^9,1250,WinLength,nOverlap,4,'linear',5,[30 70]);

if size(y,1)>size(whl,1)
    sdif = -size(whl,1)+size(y,1);
    y = y(1:end-sdif,:,:);
elseif size(y,1)<size(whl,1)
    sdif = size(whl,1)-size(y,1);
    whl = whl(1:end-sdif,:);
end

% goodi = find(whl(:,1)~=-1 & whl(:,4)=-1);
cohm = cell(3,nCh,nCh);
sm=3;
goodi = find(whl(:,1)~=-1 );

for Ch1 = 1:nCh
    for Ch2 = Ch1:nCh
       [Av, Std, Bins] = MakeAvF(whl(goodi,1:2),squeeze(y(goodi,Ch1,Ch2)),100);
%        if (Ch1==Ch2)
%            cohm{1,Ch1,Ch2} = conv2(20*log10(Av),ones(sm,sm)/sm/sm,'same');
%            cohm{2,Ch1,Ch2} = conv2(20*log10(Std),ones(sm,sm)/sm/sm,'same');
%        else
           cohm{1,Ch1,Ch2} = conv2(Av,ones(sm,sm)/sm/sm,'same');
           cohm{2,Ch1,Ch2} = conv2(Std,ones(sm,sm)/sm/sm,'same');
%        end
       cohm{3,Ch1,Ch2} = Bins;
    end
end
cnt=1;
nCh=5;
for Ch1 = 1:nCh
    for Ch2 = 1:nCh
        if (Ch2==Ch1)
            subplot(nCh+1,nCh,cnt);
            imagesc(cohm{3,Ch1,Ch2}(:,1),cohm{3,Ch1,Ch2}(:,2),log(cohm{1,Ch1,Ch2})');
            set(gca,'ydir','normal');
            mcol = median(log(cohm{1,Ch1,Ch2}(find(~isnan(cohm{1,Ch1,Ch2}(:))))));
            %caxis([-5 5]+mcol);
%             caxis([130 140]);
            axis off
            colorbar
            res = load([fn '.res.' num2str(Ch1)]);
            clu = load([fn '.clu.' num2str(Ch1)]);
            clu=clu(2:end);
            res = res(find(clu>1));
            subplot(nCh+1,nCh,cnt+nCh);
            NuPlaceField(res,whl,350,5,50,[]);
            axis off 
            colorbar
%             PlaceField(res,whl)
%             set(gca,'ydir','normal');
%             set(gca,'xdir','reverse');

        elseif (Ch2>Ch1)
            subplot(nCh+1,nCh,cnt);
            imagesc(cohm{3,Ch1,Ch2}(:,1),cohm{3,Ch1,Ch2}(:,2),(cohm{1,Ch1,Ch2})');
            set(gca,'ydir','normal');
            mcol = median(cohm{1,Ch1,Ch2}(find(~isnan(cohm{1,Ch1,Ch2}(:)))));
            %caxis([mcol-0.1 mcol+0.2]);
            axis off 
            colorbar
        end
        
            %could do pcolor instead and no reversion
            %hold on
            %contour(cohm{3,Ch1,Ch2}(:,1),cohm{3,Ch1,Ch2}(:,2),cohm{2,Ch1,Ch2}');

        cnt =cnt+1;
    end
end
