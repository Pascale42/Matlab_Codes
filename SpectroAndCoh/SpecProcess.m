%function [SpMap, Dat] = SpecProcess(Total, Index)
% auxillary function for spatial topography analisis of spectra ..
function Dat = SpecProcess(Total, varargin)

%  		Total = struct('ErpMap',ErpMap,'AvErp',AvErp,'Erps',Erps,'MaxErpCh', ...
%            MaxErpCh,'PeakAmp',PeakAmp,'CCG',Trig2Units, 'tCCG',tbin,'terp',[-nTime/3+1:nTime*2/3],	  'Latencies',Latencies*1000/SampleRate,'TrigSpec',TrigSpec,'AllSpecs',spec,'fspec',f,'tspec',t,'Channels',Channels);
[Index ] = DefaultArgs(varargin,{[1:size(Total.AllSpecs,4)]});

f = Total.fspec(:);
t = Total.tspec(:);%t =t-(t(end)-t(1))/3;
spf = find(f>5 & f<30);
spt = find(t>0.3 & t<1);
%if 0

AvSpPowMean = log(sq(mean(mean(mean(Total.AllSpecs(spt,spf,:,:),2),1),4)));
AvSpPowStd = log(sq(std(mean(mean(Total.AllSpecs(spt,spf,:,:),2),1),0,4)));
%intergrate freq. in spindle band
SpPow = log(sq(mean(Total.AllSpecs(spt,spf,:,Index),2)));
% SpPow = nT x nch x nevents
[nT nCh nEv] = size(SpPow);
%NORM SPPOW!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

%  	SpPow  = (SpPow - repmat(AvSpPowMean(:)',[nT 1 nEv]))./ ...
% 		repmat(AvSpPowStd(:)',[nT 1 nEv]);
%SpPow  = (SpPow - repmat(AvSpPowMean(:)',[nT 1 nEv]));
%SpPow = log(SpPow);
%Now find for each event and channel where max in time was
[SpPowMax SpPowLati] = max(SpPow,[],1);
SpPowMax = sq(SpPowMax);
SpPowLati = sq(SpPowLati);
[spmax spmaxch] = max(SpPowMax,[],1);
spmaxt = sq(SpPowLati(spmaxch))';
spmaxch = sq(spmaxch)';
% spmaxch - channel with max power for each event
ix = repmat(spmaxt,size(SpPow,2),1);
iy = repmat([1:size(SpPow,2)]',size(SpPow,3),1);
iz = reshape(repmat([1:size(SpPow,3)],size(SpPow,2),1),[],1);
ind = sub2ind(size(SpPow),ix,iy,iz);
MaxSpPow = reshape(SpPow(ind),size(SpPow,2),[]);

MaxSpMap = mean(MaxSpPow,2);
SpPowLats = t(spt(sq(SpPowLati(spmaxch))));
SpPowLat = median(SpPowLats);
SpMaxChs = spmaxch;

%MaxSpDiam = median(sum(MaxSpPow>repmat(std(MaxSpPow,0,1),size(MaxSpPow,1),1),1));
MaxSpDiam = PeakDiam(MaxSpPow,0.5);
%here we need to change: do some geometric center of mass weighted by spindle power

%NB!! this is not the best estimate
%	SpMaxCh = round(median(SpMaxChs));
% come back to max of the power
[dummy SpMaxCh] = max(MaxSpMap);
%keyboard

%NB!!	%SpMaxCh = round(median(sq(spmaxch))); this may not be very sorbust :(((
%SpMaxChB = bootstrp(100, 'median',SpMaxCh);
%end
AvSpPow = log(sq(mean(mean(Total.AllSpecs(spt,spf,:,Index),2),1)));
% intergral power : channel x index


SpMap = mean(AvSpPow,2);
SpMap = SpMap(:);
%keyboard
%old definition of spmaxch	%[SpMaxPow SpMaxCh] = max(SpMap);

% now the erp dat
%if isfield(Total,'Erps')
rc = 1000/1250;
nErp = size(Total.Erps,3);
for i=1:nErp
    [ErpLats(i) ErpMaxChs(i)] = max2d(sq(Total.Erps(:,:,i)));
end
ErpLats = Total.terp(ErpLats)*rc;
averp = sq(mean(Total.Erps(:,:,Index),3));
%elseif isfield(Total,'AvErp')
%	averp = Total.AvErp;
%end

goodt = find(Total.terp*rc>10 & Total.terp*rc < 500);
[ErpAmps ErpLats] = min(averp(goodt,:),[],1);
[ErpAmp ErpMaxCh] = min(sq(ErpAmps));
ErpLat = ErpLats(ErpMaxCh);
ErpMap = abs(sq(averp(goodt(ErpLat),:)));
ErpLat = Total.terp(goodt(ErpLat))*rc;

Dat = struct('SpMap',SpMap,'ErpMap',ErpMap,'MaxSpMap',MaxSpMap,'AvSpPow',AvSpPow,...
    'MaxSpPow',MaxSpPow, 'SpMaxCh',SpMaxCh,'SpMaxChs',SpMaxChs, 'SpPowLat',SpPowLat,...
    'SpPowLats',SpPowLats,'ErpAmp',abs(ErpAmp), 'ErpMaxCh',ErpMaxCh,'ErpLat',ErpLat, ...
    'AvErp',averp, 'ErpLats',ErpLats,'ErpMaxChs',ErpMaxChs,'MaxSpDiam',MaxSpDiam);
return

