	%function ccgAnalysis(DbName, Restrict, varargin)
DbName = 'mice';
%change here to get other type
Restrict = {'s','', '', '', 0 , 1};
db = RestrictDb(DbName,  Restrict);

total = GoThroughDb('Osc2Spikes',db, 0, 1, 1, 10, 30, 0,'scale',1);

%[Av, Std, Distr]=StatCell(tot,1,1,10);
fccg=[];
n = length(total);

%compute and plot means for field
for k=1:n
    fccg(:,:,:,k) = squeeze(total{k}{1}(:,1:6,1:6)); 
end
mfccg = mean(fccg,4);
errccg = std(fccg,0,4)/sqrt(n);
figure
PlotMatrix(total{1}{2}, mfccg,'b');
ForAllSubplots('hold on');
PlotMatrix(total{1}{2}, mfccg-errccg,'g');
PlotMatrix(total{1}{2}, mfccg+errccg,'g');


%now process field to units

cufccg =[];
hufccg =[];
for k=1:n
    hu = findcell(total{k}{3},'h','eq');
    cu = findcell(total{k}{3},'c','eq');
    if ~isempty(hu)
        for u=1:length(hu)
            hufccg(:,:,end+1) = squeeze(total{k}{1}(:,1:6,hu(u)));
        end
    end
    if ~isempty(cu)
        for u=1:length(cu)
            cufccg(:,:,end+1) = squeeze(total{k}{1}(:,1:6,cu(u)));
        end
    end
end
mhufccg = squeeze(mean(hufccg,3));
errhufccg = squeeze(std(hufccg,0,3)/sqrt(n));
mcufccg = squeeze(mean(cufccg,3));
errcufccg = squeeze(std(cufccg,0,3)/sqrt(n));

mufccg = cat(2, reshape(mhufccg, [], 1, 6), reshape(mcufccg, [], 1, 6));
errufccg = cat(2, reshape(errhufccg, [], 1, 6), reshape(errcufccg, [], 1, 6));
figure
PlotMatrix(total{1}{2}, mufccg,'b');
ForAllSubplots('hold on');
PlotMatrix(total{1}{2}, mufccg-errufccg,'g');
PlotMatrix(total{1}{2}, mufccg+errufccg,'g');

% hufccg = permute(hufccg, [ 4 1 2 3]);
% cufccg = permute(cufccg, [ 4 1 2 3]);
nBins = 20;
maxvalue =5;
dhuf=[]; dcuf=[]; hbins=[]; cbins=[];
 for k1 = 1:6
            be = linspace(0, min(maxvalue, max(max(hufccg(:,k1,:)))), nBins+1);
            dhuf(:,:,1,k1) = histcI(squeeze(hufccg(:,k1,:))', be);
            hbins(:,1,k1) = (be(1:end-1)+be(2:end))/2;
            
             be = linspace(0, min(maxvalue,max(max(cufccg(:,k1,:))) ), nBins+1);
            dcuf(:,:,1,k1) = histcI(squeeze(cufccg(:,k1,:))', be);
            cbins(:,1,k1) = (be(1:end-1)+be(2:end))/2;
            
%          [dhuf(:,:,1,k1) hbins(:,1,k1)] = hist(squeeze(hufccg(:,k1,:))', nBins);
%          [dcuf(:,:,1,k1) cbins(:,1,k1) ]= hist(squeeze(cufccg(:,k1,:))', nBins);

 end
 figure
%  duf = cat(3, dhuf, dcuf);
pm = SmoothMatrix(dhuf, 0.01, 0.02);
pm = reshape(pm, nBins,[], 3,2);
hbins= reshape(hbins, nBins, 3,2);
 ImageMatrix(total{1}{2},hbins, pm);
%  ForAllSubplots('ylim([0 5])');
 
 figure
 pm = SmoothMatrix(dcuf, 0.01, 0.02);
pm = reshape(pm, nBins,[], 3,2);
cbins= reshape(cbins, nBins, 3,2);
 ImageMatrix(total{1}{2},cbins, pm);
%   ForAllSubplots('ylim([0 5])');

%nmow make figure 1

figure

subplot(421)
plot(total{1}{2}, mfccg(:,1,6),'b');
hold on
plot(total{1}{2}, mfccg(:,1,6)-errccg(:,1,6),'g');
plot(total{1}{2}, mfccg(:,1,6)+errccg(:,1,6),'g');
title('delta(+) to ripples(-)');

subplot(422)
plot(total{1}{2}, mfccg(:,2,6),'b');
hold on
plot(total{1}{2}, mfccg(:,2,6)-errccg(:,2,6),'g');
plot(total{1}{2}, mfccg(:,2,6)+errccg(:,2,6),'g');
title('spindle(-) to ripples(-)');


% select some unit2filed matrixes
um=[];
um(:,:,:,1) = mufccg(:,1:2,:);
um(:,:,:,2) = mufccg(:,1:2,:)-errufccg(:,1:2,:);
um(:,:,:,3) = mufccg(:,1:2,:)+errufccg(:,1:2,:);
um = permute(um,[1 4 2 3]);

subplot(423)
plot(total{1}{2}, um(:,:,1,1));
title('delta(+) to HIP units');

subplot(424)
plot(total{1}{2}, um(:,:,1,2));
title('spindle(-) to HIP units');

subplot(425)
plot(total{1}{2}, um(:,:,2,1));
title('delta(+) to CX units');

subplot(426)
plot(total{1}{2}, um(:,:,2,2));
title('spindle(-) to CX units');

subplot(427)
plot(total{1}{2}, um(:,:,1,6));
title('ripple(-) to HIP units');

subplot(428)
plot(total{1}{2}, um(:,:,2,6));
title('ripple(-) to CX units');

ForAllSubplots('axis tight')
ForAllSubplots('grid on')

