% function CheckCCG2CluLP(FileBase, sg1, clu1, sg2, clu2, fMode, varagin)
% 
% fMode: 'c'  'd'  
% 'c_all' (need to compute AllCCGs first)
% 'd_choose' (with selected  sg1, clu1, sg2, clu2)
% 'd_all'  for all the shared bins > 25%
% 
%
% [HalfBins, BinSize] = ...
%     DefaultArgs(varargin,{21,[]});


function CheckCCG2CluLP(FileBase, sg1, clu1, sg2, clu2, fMode,varargin)

[HalfBins, BinSize] = ...
    DefaultArgs(varargin,{21,[]});

switch fMode
    
    case 'c'   % compute
        
        Par = LoadPar([FileBase '.xml']);
        BinSize = floor(Par.SampleRate/1000); % 1 ms
        
        % 1st spikegroup to load
        g1=LoadClu([FileBase '.clu.' num2str(sg1)]);
        t1=load([FileBase '.res.' num2str(sg1)]);
        
        % 2nd spikegroup to load
        g2=LoadClu([FileBase '.clu.' num2str(sg2)]);
        t2=load([FileBase '.res.' num2str(sg2)]);
        
        % cluster 1
        id1= g1 == clu1;
        t1=t1(id1);
        
        % cluster 2
        id2= g2 == clu2;
        t2=t2(id2);
        
        [ccg, t] = CCG([t1;t2],[ones(size(t1));2*ones(size(t2))], BinSize, HalfBins, ...
            Par.SampleRate,unique([ones(size(t1));2*ones(size(t2))]),'count');
        
        totccg=sum(ccg(:,1,2));
        pikccg=ccg(HalfBins+1,1,2);
        perccg=round(pikccg*100/totccg);
        
        if ~exist([FileBase '.CCG2CluLP'],'dir')
            s = sprintf('mkdir %s.CCG2CluLP', FileBase);
            system(s);
        end
        dn=[FileBase '.CCG2CluLP'];
        cd(dn)
        save([FileBase '.CCG2CluLP' num2str(sg1) '.'  num2str(clu1) '-' num2str(sg2) '.' num2str(clu2) '.mat']...
            ,'totccg','pikccg','perccg')
        cd ..
        
    case 'd'    % display
        
        dn=[FileBase '.CCG2CluLP'];
        cd(dn)
        load([FileBase '.CCG2CluLP' num2str(sg1) '.'  num2str(clu1) '-' num2str(sg2) '.' num2str(clu2) '.mat']...
            ,'totccg','pikccg','perccg','t', 'ccg')
        cd ..
        figure(7562956); bar(t,ccg(:,1,2))
        title(['Central Bin is ' num2str(pikccg) '/' num2str(totccg) ' counts for ' num2str(sg1) '.' num2str(clu1) ' vs ' num2str(sg2) ...
            '.' num2str(clu2) ' = ' num2str(perccg) '%'])
        
    case 'c_all'
        
        load([FileBase '.AllCCGs.ec.mat']);
        Par = LoadPar([FileBase '.xml']);
        BinSize = floor(Par.SampleRate/1000); % 1 ms
        
        totccg=squeeze(sum(ccg(:,:,:),1));
        pikccg=squeeze(ccg(HalfBins+1,:,:));
        perccg=round(pikccg*100./totccg);
        [rows,cols] =find(perccg >= 25);
        SameClus=[Map(rows,2:3),Map(cols,2:3),NaN(length(rows),1)];
         for n=1:length(rows)
             SameClus(n,5)=perccg(rows(n),cols(n));
         end
        
        if ~exist([FileBase '.CCG2CluLP'],'dir')
            s = sprintf('mkdir %s.CCG2CluLP', FileBase);
            system(s);
        end
        dn=[FileBase '.CCG2CluLP']; cd(dn)
        save([dn '.mat'],'totccg','pikccg','perccg','SameClus'); cd ..; clear cn
        
    case 'd_choose'
        
        load([FileBase '.AllCCGs.ec.mat'],'ccg','tbin','Map');
        Par = LoadPar([FileBase '.xml']);
        BinSize = floor(Par.SampleRate/1000); % 1 ms
        dn=[FileBase '.CCG2CluLP']; cd(dn)
        load([FileBase '.CCG2CluLP.mat'],'totccg','pikccg','perccg'); cd ..
        
        figure(420001);
        for n=1:length(sg1)
            for nn=1:length(sg2)
                [~,~, idpre] = Intersection([sg1(n) clu1(n)], Map(:,2:3));
                [~,~, idpost] = Intersection([sg2(nn) clu2(nn)], Map(:,2:3));
                
                bar(tbin,ccg(:,idpre,idpost))
                title(['Central Bin is ' num2str(pikccg(idpre,idpost)) '/' num2str(totccg(idpre,idpost)) ' counts for ' ...
                    num2str(sg1(n)) '.' num2str(clu1(n)) ' vs ' num2str(sg2(nn)) '.' num2str(clu2(nn)) ' = ' num2str(perccg(idpre,idpost)) '%'])
                waitforbuttonpress
            end
        end
        
    case 'd_all'
        
        load([FileBase '.AllCCGs.ec.mat']);
        Par = LoadPar([FileBase '.xml']);
        BinSize = floor(Par.SampleRate/1000); % 1 ms
        dn=[FileBase '.CCG2CluLP']; cd(dn)
        load([dn '.mat'],'totccg','pikccg','perccg'); cd ..; clear dn
        [rows,cols] =find(perccg >= 25);
        figure(420002);
        for n=1:length(rows)
            subplot 131; bar(tbin,ccg(:,rows(n),rows(n))); axis tight; title(num2str(Map(rows(n),2:3)))
            subplot 132; bar(tbin,ccg(:,cols(n),cols(n))); axis tight;title(num2str(Map(cols(n),2:3)))
            subplot 133;  bar(tbin,ccg(:,rows(n),cols(n))); axis tight;
            title(['Central Bin is ' num2str(pikccg(rows(n),cols(n))) '/' num2str(totccg(rows(n),cols(n))) ' counts for ' ...
                num2str(Map(rows(n),2:3)) ' vs ' num2str(Map(cols(n),2:3)) ' = ' num2str(perccg(rows(n),cols(n))) '%'])
            waitforbuttonpress
        end
        
end


