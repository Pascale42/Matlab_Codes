function out = WithinRangesLong(x, Ranges, varargin)
%out = WithinRangesLong(x, Ranges, <RangeLabel>, <Mode>, <BlockSize>, <SparseOut>)
%
%It does the same as WithinRanges(x, Ranges, RangeLabel, Mode)
%but works much faster for high (>100) number of Ranges, especially when RangeLabel is provided.
%It just runs WithinRangesLong in a loop over smaller blocks of Ranges of length BlockSize (default = 50).
%Note: RangeLabel must be a vector 1:nRanges!
%
% Evgeny Resnik
% version 12.09.2012
% added sparse output as option to save memory. 17.12.2013 AS
%      out = WithinRangesLong(spk.Tlfp(GoodSpk), BurstTimeSpan(GoodBursts,:) , [1:length(GoodBursts)]');

%NOTE: when RangeLabel=[] - is not implemented yet.


nRanges = size(Ranges,1);

% Parse input parameters
[ RangeLabel, Mode,  BlockSize, SparseOut ] = DefaultArgs(varargin,{ ones(nRanges,1), 'matrix', 50, 0 });

%Check whether non-empty RangeLabel is provided
if isequal(RangeLabel, ones(nRanges,1) )
    %error('Case when RangeLabel=[] is not implemented yet!')
    IfRangeLabel=0;
else
    IfRangeLabel=1;
end

nBlocks =  fix(nRanges/BlockSize);
RestSize = nRanges - BlockSize*nBlocks;

switch IfRangeLabel
    case 0   %RangeLabel is empty or not provided       
     
%         for b=1:nBlocks
%             BlockInd = [ (b-1)*BlockSize+1 : b*BlockSize ];
%             out0(:,b) = WithinRanges(x, Ranges(BlockInd,:) );
%         end        
%         %rest
%         BlockInd = nRanges-RestSize+1 : nRanges;
%         out0(:, end+1) = WithinRanges(x, Ranges(BlockInd,:) );
%         out = double(sum(out0,2)>0);
        
        %just do it as the original function because it is fast enough
        out = WithinRanges(x, Ranges);
        
        
    case 1  %non-empty RangeLabel is provided        
       
        
        %out = []; slow version with high length(x) and nRanges
        if ~SparseOut
            out = NaN*ones(length(x), nRanges);
        else
            maxnz = 10*length(x); 
            out=  spalloc(length(x), nRanges,maxnz);
        end
        
        for b=1:nBlocks
            %disp(['block-' num2str(b)])
            BlockInd = [ (b-1)*BlockSize+1 : b*BlockSize ];
            out0 = WithinRanges(x, Ranges(BlockInd,:) , [1:length(BlockInd)]', Mode);
            out(:,BlockInd) = out0;
            
            %out = cat(2,out, out0); %slow version with high length(x) and nRanges
        end
        
        %rest
        if RestSize>0
            BlockInd = nRanges-RestSize+1 : nRanges;
            out0 = WithinRanges(x, Ranges(BlockInd,:) , [1:length(BlockInd)]', Mode);
            out(:,BlockInd) = out0;
            %out = cat(2,out, out0); %slow version with high length(x) and nRanges
        end
        
end %switch IfRangeLabel



