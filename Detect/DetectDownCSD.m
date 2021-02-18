%detection of down state in CSD 
%function out = DetectDownCSD(FileBase, fMode, Channels,FreqRange)
% fMode = 'compute'
function out = DetectDownCSD(FileBase, varargin)
%FileBase = 'sm9603_228-297';
Par = LoadPar([FileBase '.xml']);
eFs= 1250;
[fMode, Channels, FreqRange] = DefaultArgs(varargin,{'compute',load([FileBase '.ecsink']),[100 600]});
State ='SWS';
[RepChan Tetrode] = RepresentChan(FileBase);
Period  = load([FileBase '.sts.' State]);
nT = FileLength([FileBase '.eeg'])/Par.nChannels/2;

csd = LoadBinary([FileBase '.csd'], Channels, Par.nChannels)';
scsd = ButFilter(csd,4,10/625,'low');
fcsd = ButFilter(csd,4,[50 130]/625,'bandpass');
fcsd = Filter0(ones(100,1)/100,abs(fcsd));
fcsd = unity(fcsd);
nT = length(fcsd);

uthr =-1.5;
dthr = -1.5;
[ucr dcr] = SchmittTrigger(fcsd,uthr,dthr);

%keep only within periods
[ucr gi] = SelectPeriods(ucr,Period,'d',1);
dcr = dcr(gi);
[dcr gi] = SelectPeriods(dcr,Period,'d',1);
ucr = ucr(gi);

%filter crossings to remove short gaps
duint = [Inf; dcr(2:end-1)-ucr(1:end-2)];
mergei = find(duint<1.25*20);
ucr(mergei) = ucr(mergei+1);
ucr(mergei+1)=[];
dcr(mergei+1)=[];

len= (ucr-dcr)./1.25;

if 1
    gi = find(len>50 & len<1000);
    len = len(gi); ucr=ucr(gi); dcr = dcr(gi);
end

lms = LocalMinima(scsd);
[lms lmsi ucri] = NearestNeighbour(lms,ucr,'right');
lms = unique(lms);
ncr = length(lms);

Win = [625 625];
Winlen = sum(Win)+1;
[ssegc compl] = GetSegs(scsd,lms-Win(1),Winlen,[]);
lms = lms(compl); ncr = length(lms);

seg = GetSegs(fcsd,lms-Win(1),Winlen,[]);
fseg = ButFilter(seg,4,10/625,'low');

%[dummy maxi] = max(ssegc(1:Win(1),:));
[dummy maxi] = max(fseg(1:Win(1),:));

UDU = [maxi(:)-Win(1)+lms(:) lms(:)];

%now filter by length
len = diff(UDU,1,2)/1.25;
gi = find(len>50 & len<Win(1));
UDU=UDU(gi,:); lms = lms(gi); ncr =length(lms);len = len(gi);
if 0
uduind = WithinRanges([1:nT],UDU,[1:ncr],'vector');
gi = uduind>0;
mgam_down = Accumulate(uduind(gi)',fcsd(gi),ncr);
vgam_down = Accumulate(uduind(gi)',fcsd(gi).^2,ncr);
n_down = Accumulate(uduind(gi)',1,ncr);
mgam_down = mgam_down./n_down;
vgam_down = vgam_down./n_down-mgam_down.^2;
% 
% seg = GetSegs(fcsd,ucr-625,1251,[]);
% segc = GetSegs(csd,ucr-625,1251,[]);
% lns = [dcr-ucr ucr-ucr]+625;
end
segc = GetSegs(csd,lms-Win(1),Winlen,[]);
seg = GetSegs(fcsd,lms-Win(1),Winlen,[]);
PlotTraces(segc,[-Win(1):Win(2)]/1.25,1250,100);

%useg =abs(unity(seg(1:625,:)));

%fsegc = ButFilter(segc,4,10/625,'low');
% if 0
% for i=1:length(lns)
% %    lm = LocalMinima(-fseg(1:625,i));
%      mingam = min(fseg(1:625,i));
%      maxgam =fseg(625,i);
%      b= (mingam+maxgam)/2;
%      a=[];
%      while isempty(a) 
%          a=SchmittTrigger(flipud(fseg(1:625,i)),b,b);
%          b=(mingam+b)/2;
%      end
%      downbeg(i) = max(625-a);
% end
% end


keyboard
[Res Clu ] = ReadEl4CCG(FileBase);
Res = round(Res/16);
[Res ind] = SelectPeriods(Res,Period,'d',1);
Clu=Clu(ind);
CluLoc = load([FileBase '.cluloc']);
dg = find(CluLoc(:,4)==3);
ca1 = find(CluLoc(:,4)==1);
ca3 = find(CluLoc(:,4)==2);

% imagesc(fsegc(:,si)')
% segind = WithinRanges(Res(:),[dcr ucr],[1:ncr],'vector');
% 
% gspk = find(segind>0);
% resinseg = Res(gspk)-dcr(si(segind(gspk)'));
% cluinseg = Clu(gspk);
% segind = segind(gspk);

[T,G] = CatTrains({lms,Res},{1,Clu});
[T,G1] = CatTrains({lms,Res},{[1:ncr]',Clu});
[ccg tbin pairs] = CCG(T,G,Winlen,0,1250,[1:max(G)],'count');
MyPairs = pairs(G(pairs(:,1))==1 & G(pairs(:,2))>1,:);
SegRes = diff(T(MyPairs),1,2);
SegClu = G(MyPairs(:,2));
SegInd = G1(MyPairs(:,1));

DentInd = find(ismember(SegClu,dg));
DentDown = (SegRes(DentInd) > maxi(SegInd(DentInd))'-625) & (SegRes(DentInd) <0);
DentCnt = Accumulate(SegInd(DentInd),double(DentDown),ncr);
NoDent = find(DentCnt==0);
SilentI = ismember(SegInd,NoDent);
SegIndS = SegInd(SilentI);
SegCluS = SegClu(SilentI);
SegResS = SegRes(SilentI);

%NoDent = [1:ncr];

dgres= Res(ismember(Clu,dg));
ca3res= Res(ismember(Clu,ca3));
ca1res= Res(ismember(Clu,ca1));
GoodSeg = intersect(NoDent,si(1:200));
[T1,G1] = CatTrains({lms(GoodSeg),lms(GoodSeg)-625+maxi(GoodSeg)',ca1res,ca3res,dgres},{1,2,3,4,5});
figure(323232);clf
CCG(T1,G1,20,40,1250,[1:max(G1)],'count');
ForAllSubplots('xlim([-500 500]);');
ForAllSubplots('Lines(0,[],''r'');');

[s si] = sort(len);
[ss ssi] = sort(si);

ca1i = ismember(SegCluS,ca1);
ca3i = ismember(SegCluS,ca3);
dgi = ismember(SegCluS,dg);
figure(3232);clf
plot(SegResS(ca1i)/1.25,ssi(SegIndS(ca1i)),'.b','MarkerSize',4);
hold on
plot(SegResS(ca3i)/1.25,ssi(SegIndS(ca3i)),'.g','MarkerSize',6);
plot(SegResS(dgi)/1.25,ssi(SegIndS(dgi)),'.r','MarkerSize',6);


 figure(365656);clf
 [dummy dymii Ca3Clu] = unique(SegCluS(ca3i));
 scatter(SegResS(ca3i)/1.25,ssi(SegIndS(ca3i)),15,Ca3Clu,'filled');
 hold on
 plot(SegResS(dgi)/1.25,ssi(SegIndS(dgi)),'.k','MarkerSize',4);

 
[T,G] = CatTrains({lms(GoodSeg),lms(GoodSeg)-625+maxi(GoodSeg)',Res},{1,2,Clu}); 
[T,G] = CatTrains({nlms,Res},{1,Clu}); 
[ccg tbin] = CCG(T,G,20,40,1250,[1:max(G)],'count');
figure
for i=1:max(G)
    subplotfit(i,max(G));
    bar(tbin,ccg(:,1,1+i));axis tight
    title(num2str(CluLoc(i,[1 2 4])));
end
ForAllSubplots('xlim([-400 400]);');
ForAllSubplots('Lines(0,[],''r'');');

% for i=1:length(lns)
%     clf
%     plot([unity(fsegc(:,i)) seg(:,i)]);axis tight
%     hold on
%     Lines(lns(i,:),[],'r');
%     waitforbuttonpress;
% end
% 
% 
