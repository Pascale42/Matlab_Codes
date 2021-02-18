function MonoModExample(FileBase,varargin)
if ~FileExists([FileBase '.mono-c'])
    return
end
load([FileBase '.mono-c'],'-MAT');
if isempty(mono.From.ElClu)
    return
end
% RootDir = pwd;
% cd(FileBase);
Par = LoadPar([FileBase '.xml']);
El = find(strcmp(Par.ElecLoc,'c'));

phmod = load([FileBase '.PhaseFrAmpAnal.mat']);
phmodElClu = CatStruct(phmod.OutArgs,{'El','Clu'});

load([FileBase '.s2s'],'-MAT');
load([FileBase '.NeuronQuality.mat']);
OutArgs = CatStruct(OutArgs);
mycnq = find(ismember(OutArgs.ElNum,El));
nq = SubsetStruct(OutArgs, mycnq,1);

[el clu spkwdths ie] = textread([FileBase '.type-c'],'%d %d %s %d');

n=size(mono.From.ElClu,1);
mutind = [];
muted = zeros(n,1);
for ii=1:n
    if mono.To.Type==0 | muted(ii)==1
        continue;
    else
        tmp = find(mono.From.ElClu(:,1)==mono.To.ElClu(ii,1) & mono.From.ElClu(:,2)==mono.To.ElClu(ii,2) ...
            & mono.To.ElClu(:,1)==mono.From.ElClu(ii,1) & mono.To.ElClu(:,2)==mono.From.ElClu(ii,2));
        if ~isempty(tmp)
            mutind(end+1) = tmp;
            muted(tmp) = 1;
        end
    end
end
mono.From.ElClu(mutind,:)=[];
mono.To.ElClu(mutind,:)=[];
mono.From.Type(mutind)=[];
mono.To.Type(mutind)=[];

n=size(mono.From.ElClu,1);
for ii=1:n
    %now get index from phmod
    fromi = find(phmodElClu.El==mono.From.ElClu(ii,1) & phmodElClu.Clu==mono.From.ElClu(ii,2));
    toi = find(phmodElClu.El==mono.To.ElClu(ii,1) & phmodElClu.Clu==mono.To.ElClu(ii,2));
    myphmod = phmod.OutArgs([fromi toi]);
    if isempty(toi) | isempty(fromi)
        continue;
    end
    figure(9843);clf
    subplot(3,3,1)
    bar([myphmod(1).phbin myphmod(1).phbin+360],[myphmod(1).phhist myphmod(1).phhist]);
    hold on
    stairs([myphmod(1).phbin myphmod(1).phbin+360],2*Filter0(gausswin(3)/3,[myphmod(1).phhist myphmod(1).phhist]),'r');
    axis tight
    title([num2str(mono.From.ElClu(ii,:)) ', Type ' num2str(mono.From.Type(ii))]);    
    
    subplot(3,3,4)
    bar([myphmod(2).phbin myphmod(2).phbin+360],[myphmod(2).phhist myphmod(2).phhist]);
    hold on
    stairs([myphmod(2).phbin myphmod(2).phbin+360],2*Filter0(gausswin(3)/3,[myphmod(2).phhist myphmod(2).phhist]),'r');
    axis tight
    title([num2str(mono.To.ElClu(ii,:)) ', Type ' num2str(mono.To.Type(ii))]);
    
    fr1 = myphmod(1).n/myphmod(1).StateTime;
    fr2 = myphmod(2).n/myphmod(2).StateTime;
   
    %now the ccgs
    fromi = find(s2s.ElClu(:,1)==mono.From.ElClu(ii,1) & s2s.ElClu(:,2)==mono.From.ElClu(ii,2));
    toi = find(s2s.ElClu(:,1)==mono.To.ElClu(ii,1) & s2s.ElClu(:,2)==mono.To.ElClu(ii,2));
    
    subplot(3,3,2)    
    bar(s2s.tbin,sq(s2s.ccg(:,fromi,fromi,4))'); xlim([-50 50]);
    title(['FirRate = ' sprintf('%2.2f',fr1) ]);
    subplot(3,3,5)    
    bar(s2s.tbin,sq(s2s.ccg(:,toi,toi,4))'); xlim([-50 50]);
    title(['FirRate = ' sprintf('%2.2f',fr2) ]);
    
    statei = find(strcmp(s2s.State,myphmod(1).State));
    if fromi<toi
        ccg = sq(s2s.ccg(:,fromi,toi,4));
        ccgs = sq(s2s.ccg(:,fromi,toi,statei));
    else
        ccg = flipud(sq(s2s.ccg(:,toi,fromi,4)));
        ccgs = flipud(sq(s2s.ccg(:,toi,fromi,statei)));
    end
    subplot(3,3,7)    
    bar(s2s.tbin,ccg); xlim([-30 30]);
    
    subplot(3,3,8)    
    bar(s2s.tbin,Filter0(ones(5,1)/5,ccgs)); xlim([-250 250]); %axis tight
   
    %now the waveshape
    fromi = find(nq.ElNum==mono.From.ElClu(ii,1) & nq.Clus==mono.From.ElClu(ii,2));
    toi = find(nq.ElNum==mono.To.ElClu(ii,1) & nq.Clus==mono.To.ElClu(ii,2));
    
    subplot(3,3,3)
    plot([-32:31]./40,nq.AvSpk(fromi,:)); axis tight
    title(['IsPos = ' num2str(nq.IsPositive(fromi)) ', ' spkwdths{fromi}]);
    subplot(3,3,6)
    plot([-32:31]./40,nq.AvSpk(toi,:)); axis tight
    title(['IsPos = ' num2str(nq.IsPositive(toi)) ', ' spkwdths{toi}]);    
    suptitle([FileBase ' ' myphmod(1).State]);
    
    infostr = sprintf('%s , PRE: pval=%2.3f, th0=%2.1f, POST: pval=%2.3f, th0=%2.1f', ...
        FileBase, myphmod(1).pR, myphmod(1).th0*180/pi,myphmod(2).pR, myphmod(2).th0*180/pi);
  [ dummy dymmmi b] = PointInput(1);
  switch b
      case 1
          continue;
      case 2
          return
      case 3
          reportfig(gcf,'MonoModExamples',0,infostr);
  end
      
end


%cd(RoodDir);
