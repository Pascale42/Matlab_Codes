
FileBase='';

%% Clus to pairs to Pairs

AllCCGs(FileBase, 'hpc', 'c'); % compute

AllCCGs(FileBase, 'ec', 'd1')  % display

AllCCGs(FileBase, 'old', 'd2',[x x]); % for one Cluster

AllCCGs(FileBase, 'old', 'd3', CluPre(x:x,:));

% pairs created manually after visual detection of suspicious pairs after AllCCGs.m

[dum,dum,Map]=LoadCluRes(FileBase,[1:4]); clear dum  % if not computed

[dum,dum,Map]=load([FileBase 'CluRes.mat']); clear dum

Pairs=zeros(size(pairs,1),4); 

for n=1:size(pairs,1)
[dum,dum, idclu] = intersect(pairs(n,1), Map(:,1));
Pairs(n,[1 2])= Map(idclu,[2 3]);

[dum,dum, idclu] = intersect(pairs(n,2), Map(:,1));
Pairs(n,[3 4])= Map(idclu,[2 3]);
clear dum
end
clear idclu n

save([FileBase '.CCGJ.GoodPairs.mat'],'Pairs', 'Map')

%% CCG Jitter

tic
for n=1:size(Pairs, 1)
    CCG_Jitter(FileBase,Pairs(n,1:2),Pairs(n,3:4),'c');
end
toc


for n=1:size(Pairs, 1)
    CCG_Jitter(FileBase,Pairs(n,1:2),Pairs(n,3:4),'d');
    xlabel(num2str(n))
    waitforbuttonpress
end


save([FileBase '.CCGJ.GoodPairs.mat'],'Pairs', 'Map', 'CnxCt')

