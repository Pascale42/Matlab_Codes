% [Out xr yr Count] = CCGEvolutionOverlap(Res, Clu, BinSize, HalfBins, SampleRate, GSubset, ChunkSize, Step);
% this one unlike the CCGEvolution doesnot take Segments and just computes
% CCGs in overlaping windows
% divides time into chunks of length ChunkSize (in samples) and then
% computes a CCG in each one, and displays them as a 2d plot
%
% optional outputs Out is matrix of outputs Out(CCGBin, Chunk, Cell1, Cell2) for Cell2>Cell1
% xr and yr are for if you want to display it: imagesc(xr, yr, Out(:,:,Cell1, Cell2));

function [Out, xr, yr, OutCnt]  = CCGEvolutionOverlap(T, G, varargin)

Eps = 1;

[BinSize, HalfBins, SampleRate, GSubset, ChunkSize,Step] = DefaultArgs(varargin,...
    {200,  100,     20000,      [1:max(G)], 20000,   10000});

% compute CCG
[out, t, Pairs] = CCG(T, G, BinSize, HalfBins, SampleRate, GSubset);

%find how many chunks we have to skip to get next non-overlapping -
% this will define the vectorized computation scheme: we do the same as for
% nonoverlapped for several (how many see below) series of chunks. plus
% some extra at the end
maxT =max(T);
minT =min(T);

for i1=1:length(GSubset)
    Cell1 = GSubset(i1);
    for i2=i1+1:length(GSubset) % excluded autocorrelogram for now
        Cell2=GSubset(i2);
        %fprintf('Doing %d %d\n', Cell1, Cell2);

		MyPairs = Pairs(G(Pairs(:,1))==Cell1 & G(Pairs(:,2))==Cell2, :);
		nPairs = size(Pairs,1);
        
		SpkInd = [];
        %get index of first and second elements of MyPairs from FitOverlpedBins
        [ChunkSegs SpkInd1 SegInd1] = FitOverlapedBins(ChunkSize,Step, T(MyPairs(:,1)), minT);
        [ChunkSegs SpkInd2 SegInd2] = FitOverlapedBins(ChunkSize,Step, T(MyPairs(:,2)), minT);
        nChunks = size(ChunkSegs,1);
        
        s1 = sparse(SpkInd1,SegInd1,1,nPairs,nChunks);  
        s2 = sparse(SpkInd2,SegInd2,1,nPairs,nChunks); 
        s = and(s1,s2);   
         
        %now get only those pairs that fall in the same segment
        [GoodPairs GoodSegs] = find(s);
                
        %here is in what chunk they fall
        %WhichChunk = SegInd1(GoodPairs);
        %WhichGrp = SegGrpInd(PairInd1(GoodPairs));
        
        % MATLAB SUCKS!!!! THIS WHOLE INDEXING THING IS A PILE OF SHIT
        % this finds the index of the CCg bin where the pairs fall
        if size(MyPairs,1)>1
    		DiffBin = round(diff(T(MyPairs(GoodPairs,:)),1,2)/BinSize) + HalfBins + 1;
        else
            DiffBin = round(diff(T(MyPairs(GoodPairs)))/BinSize) + HalfBins + 1;
        end
        
		if any(DiffBin<=0 | DiffBin>1+2*HalfBins)
            warning('Time difference Out of range')
            DiffBin(DiffBin<=0) = 1;
            DiffBin(DiffBin>1+2*HalfBins) = 1+2*HalfBins;
        end
        
        %get binary matrix of first and seecond spike in a pair to
        %belong to the segment

        if size(MyPairs, 1) > 0
            Count = Accumulate([DiffBin, GoodSegs], 1, [2*HalfBins+1, nChunks ]);
        else
            Count = zeros(2*HalfBins+1, nChunks);
        end

        % scale it.  How many coincidences are expected?
        
        n1 = full(sum(s1,1));
        n2 = full(sum(s2,1));
        %n1 = Accumulate(SegInd(ismember(SpkInd,find(G==Cell1))),1,nChunks);
        %n2 = Accumulate(SegInd(ismember(SpkInd,find(G==Cell2))),1,nChunks);
        
        Expected = n1.*n2/(ChunkSize/BinSize);
        gi = find(Expected>0);
        Scales = zeros(2*HalfBins+1,nChunks);
        Scaled(:,gi) = (Count(:,gi)+Eps) ./(repmat(Expected(gi),2*HalfBins+1,1)+Eps);
        
        %^        Scaled(:,zeroi) = NaN;

        if nargout>=1
            Out(:,:,i1,i2) = Scaled;
            Out(:,:,i2,i1) = Scaled;
                                 
            xr = (ChunkSegs(:,1)+ChunkSize/2)/SampleRate;
            yr = (-HalfBins:HalfBins)*BinSize/SampleRate*1000;
            
            if nargout>3
                OutCnt(:,:,i1,i2) = Count;
                OutCnt(:,:,i2,i1) = Count;
            end
            
        else
			cla; hold off
			imagesc((0:nChunks-1)*Step/SampleRate, (-HalfBins:HalfBins)*BinSize/SampleRate*1000, log10(Scaled));
			%caxis([-1 1])
			hold on; plot(xlim, [25 25], 'k:', xlim, [-25 -25], 'k:');
			colorbar
			xlabel('Time (s)');
			ylabel('dt (ms)');
            drawnow
        end
    end
end
return




if 0

    
%TO BE CONTINUED IN THE GOOD MOOD FOR VECTORIZATION! :))

%ChunkSegs = OverlapSegments(minT, maxT, ChunkSize, Step,0);
%[ChunkSegs SpkInd SegInd SegGrpInd] = FitOverlapedBins(ChunkSize,Step, T,
%minT);^% actually we don't need the group indexing now - vectorization
%class!!
%[ChunkSegs SpkInd SegInd] = FitOverlapedBins(ChunkSize,Step, T, minT);


%nChkGps = floor(ChunkSize/Step)+1;

    
    %DUMB WAY - FOR LOOP!!!
Segs = OverlapSegments(minT, maxT, BinSize, Step, 0);
nChunks = size(Segs,1);

for i1=1:length(GSubset)
    Cell1 = GSubset(i1);
    for i2=i1:length(GSubset)
        Cell2=GSubset(i2);
        %fprintf('Doing %d %d\n', Cell1, Cell2);

        MyPairs = Pairs(G(Pairs(:,1))==Cell1 & G(Pairs(:,2))==Cell2, :);

        for c=1:nChunks
            myi1 =T(MyPairs(:,1))>Segs(c,1) & T(MyPairs(:,1))<Segs(c,2);
            myi2= T(MyPairs(:,2))>Segs(c,1) & T(MyPairs(:,2))<Segs(c,2);
            myi = find(myi1 & myi2);
            % gives indexes of Pairs that both fall in chunk #c

            if isempty(myi)
                Scaled(:,c) = zeros(2*HalfBins+1,1);
            else

                if size(MyPairs,1)>1
                    DiffBin = round(diff(T(MyPairs(myi,:)),1,2)/BinSize) + HalfBins + 1;
                else
                    DiffBin = round(diff(T(MyPairs(myi,:)))/BinSize) + HalfBins + 1;
                end

                if any(DiffBin<=0 | DiffBin>1+2*HalfBins)
                    warning('Time difference Out of range')
                    DiffBin(DiffBin<=0) = 1;
                    DiffBin(DiffBin>1+2*HalfBins) = 1+2*HalfBins;
                end

                Count = Accumulate(DiffBin,1,2*HalfBins+1);
                
                n1 = sum(myi1); n2 = sum(myi2);
                
                Expected = n1.*n2/(ChunkSize/BinSize);

                Scaled(:,c) = (Count+Eps) ./(repmat(Expected',2*HalfBins+1,1)+Eps);
                if nargout>3
                    TotCnt(:,c) = Count;
                end
            end

        end
         Out(:,:,i1,i2) = Scaled;
         Out(:,:,i2,i1) = Scaled;
         if nargout>3
            OutCnt(:,:,i1,i2) = TotCount;
            OutCnt(:,:,i2,i1) = TotCount;
         end
         xr = (0:nChunks-1)*ChunkSize/SampleRate;
         yr = (-HalfBins:HalfBins)*BinSize/SampleRate*1000;
        
    end
end




%return

 end








% test code
% first 100s of 10hz independent
res1 = rand(2000, 1)*20000*100;
clu1 = [ones(1000,1);2*ones(1000,1)];

% then 100s of 10hz correlated with 25ms jitter
res2 = (1+rand(1000,1))*20000*100; res2 = [res2;res2+randn(1000,1)*25*20];
clu2 = clu1;

CCGEvolutionOverlap([res1;res2], [clu1;clu2], 80, 50, 20000, [1 2], 20000*20,20000*5);


% old way of computing using WithinRanges - very memory consuming for large
% number of spikes and bins
% for chg=1:nChkGps
%             % which chunk does each spike belong to?
%             myChkInd = chg:nChkGps:nChunks;
%             nmyChks =length(myChkInd);
%             myChunks = ChunkSegs(myChkInd,:);
%             %get binary matrix of first and seecond spike in a pair to
%             %belong to the segment
%             IfInChunk1= WithinRanges(T(MyPairs(:,1)), myChunks,[1:nmyChks],'matrix');
%             IfInChunk2= WithinRanges(T(MyPairs(:,2)), myChunks,[1:nmyChks],'matrix');
%             IfInChunk = and(IfInChunk1,IfInChunk2); % want both to belong
%             [GoodPairs SpkChunk] = find(IfInChunk);
%             %SpkChunk = myChkInd(
% 
%             if size(MyPairs, 1) > 0
%                 Count = Accumulate([DiffBin(GoodPairs), SpkChunk], 1, [2*HalfBins+1, nmyChunks]);
%             else
%                 Count = zeros(2*HalfBins+1, nmyChunks);
%             end
% 
%             % scale it.  How many coincidences are expected?
%             %n1 = Accumulate(SpkChunk(G(GoodPairs)==Cell1 & SpkChunk>0),1,nmyChunks);
%             %n2 = Accumulate(SpkChunk(G(GoodPairs)==Cell2 & SpkChunk>0),1,nmyChunks);
%             n1 = sum(IfInChunk1);
%             n2 = sum(IfInChunk2);
%             Expected = n1.*n2/(ChunkSize/BinSize);
% 
%             Scaled(:, myChkInd) = (Count+Eps) ./(repmat(Expected,1,2*HalfBins+1)+Eps);
%         