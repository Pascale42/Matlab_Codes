function [out, t, Pairs] = CCG_part(T, G, BinSize, HalfBins, SampleRate, GSubset, Normalization, Epochs)
% constructs multiple cross and Auto correlogram for subset of pairs
% usage: [ccg, t, pairs] = CCG_part(T, G, BinSize, HalfBins, SampleRate, GSubset, Normalization, Epochs)
%
% T gives the time of the events, (events need not be sorted by TIME)
% G says which one is in which group
% BinSize gives the size of a bin in input units (i.e. not scaled by SampleRate)
% HalfBins gives the number of bins on each side of 0 - so the total is 1+2*HalfBins
% SampleRate is for x-axis scaling only.  It defaults to 20000
% GSubset says which groups to plot the CCGS of (defaults to all but group 1)
% If GSubset is m x 2 matrix then it contains list of pairs btw which to
% compute CCG, this will speed up the computation when you only need some
% pairs and don't want to loop by yourself. GSubset can be a square binary
% matrix, with 1s at pairs to be computed, zeros - otherwise.
% Normalization indicates the type of y-axis normalization to be used.  
% 'count' indicates that the y axis should show the raw spike count in each bin.
% 'hz' will normalize to give the conditional intensity of cell 2 given that cell 1 fired a spike (default)
% 'hz2' will give the joint intensity, measured in hz^2.
% 'scale' will scale by both firing rates so the asymptotic value is 1.  This gives you
%   the ratio of the coincidence rate to that expected for uncorrelated spike trains
%
% optional input Epochs allows you to compute from only spikes in certain epochs
% and it will bias-correct so you don't see any triangular shape. 
% Warning: if gaps between epochs are shorter than total CCG length, this will mess
% up the edges.
%
% The output array will be 3d with the first dim being time lag and the second two 
% specifying the 2 cells in question (in the order of GSubset)
% If there is no output specified, it will plot the CCGs
%
% This file calls a C program so your CCG is computed fast.
% to use it, you need to compile mex file from CCGHeart.c
% run : mex -v CCGHeart_part.c -o CCGHeart_part_71 (if you are using matlab7.1
% ...  _72 for matlab7.2, etc .. you understand)
% different architectures will produce different extensions of mex file,
% also different versions of matlab link mex file to different libraries
% they are mostly taken into account in the code of CCG.m, but if your
% version or archicture is different from those - modify CCGFun string to
% match the name of the mex file your compiler generated.
% optional output t gives time axis for the bins in ms
% optional output argument pairs gives a nx2 array with the indices of the spikes
% in each train that fall in the CCG.
% originally written by Ken Harris 
% small aditions Anton Sirota
if nargin<5
	SampleRate = 20000;
end
if nargin<6
	GSubset = unique(G);
%	GSubset = setdiff(GSubset, 1);
end
if min(size(GSubset))==2 & any(GSubset(:)>1)
    Pairs2Use = GSubset;
    GSubset = unique(GSubset(:));
end

if size(GSubset,1) == size(GSubset,2) & all(ismember(GSubset(:),[0 1]))
    %then GSubset is a binary square matrix 
    [Pairs2Use(:,1) Pairs2Use(:,2) v] = find(logical(GSubset));
    GSubset = unique(Pairs2Use(:));
end

if nargin<7
	Normalization = 'hz';
end;
if nargin<8
    Epochs = [];
end

if length(G)==1
	G = ones(length(T), 1);
	GSubset = 1;
	nGroups = 1;
else
	nGroups = length(GSubset);
end;


% Prepare Res and Clu arrays.
G=G(:);
T=T(:);

if ~isempty(Epochs)
    Included = find(ismember(G,GSubset) & isfinite(T) & WithinRanges(T,Epochs));
    
    % check gaps between epochs are not too short
    GapLen = Epochs(2:end,1) - Epochs(1:(size(Epochs,1)-1),2);
    TooShort = find(GapLen<BinSize*(HalfBins+.5));
    if ~isempty(TooShort)
        fprintf('WARNING: Epochs ');
        fprintf('%d ', TooShort);
        fprintf('are followed by too-short gaps.\n');
    end
 

else
    Included = find(ismember(G,GSubset) & isfinite(T));
    Epochs = [];%[min(T)-1 max(T)+1];
  
end
Res = T(Included);
% if no spikes, return nothing
if length(Res)<=1
    nBins = 1+2*HalfBins;
    out = zeros(nBins, nGroups, nGroups);
    t = 1000*(-HalfBins:HalfBins)*BinSize/SampleRate;
    Pairs = [];
    return
end

% To make the Clu array we need an indexing array, which SUCKS!
G2Clu = full(sparse(GSubset,1,1:nGroups));
Clu = G2Clu(G(Included));
nPresent = max(Clu);

if exist('Pairs2Use','var')
    %now create the pair selection matrix:
    Pairs2Use = G2Clu(Pairs2Use); % renumber as well
    bUsePairs = zeros(nPresent,nPresent);
    bUsePairs(Pairs2Use(:,1),Pairs2Use(:,2)) = 1;
end

% sort by time 
[Res ind] = sort(Res);
Clu = Clu(ind);

% Now call the C program...
%check for version and call appropriatly compiled mex function 
% (libraries differ, so need to check)
ver = version; 
[dd host ] = system('hostname');
%new smart way: compile on each platform and give name:
% first 3 symbols of the version and mex extension, conveniently
% returned by matlab e.g. CCGHeart_7-1.mexext
verstr = ver([1 3]);
if exist('Pairs2Use','var')
    CCGFun = ['CCGHeart_part_' verstr ];
end
if strfind(host,'urethane') & str2num(ver(1:3))==7.2
    % god knows why urethane does not compile normally and crashes
    error('god knows why (some libraries bug?) on urethane CCGHeart compiled under matlab 7.2 does crashes, try 7.1 or 7.3');
    %         CCGFun = ['CCGHeart_urethane'];
end
%fprintf('using function name %s for version=%d mexext=%s\n',CCGFun,ver,mexext);

% call the program
nSpikes = length(Res);
if nargout>=3
    % fixing the bug of CCGHeart when no spikes fall withing HalfBins (even
    % for autocorrelogram
    if min(diff(Res))<=BinSize*HalfBins
        if exist('Pairs2Use','var')
            [Counts RawPairs] = feval(CCGFun,Res, uint32(Clu), BinSize, uint32(HalfBins),uint32(bUsePairs));
            rsRawPairs = reshape(RawPairs, [2 length(RawPairs)/2])';
        else
            [Counts RawPairs] = feval(CCGFun,Res, uint32(Clu), BinSize, uint32(HalfBins));
            rsRawPairs = reshape(RawPairs, [2 length(RawPairs)/2])';
        end
    else
        warning('pairs cannot be computed - no overlap between spikes withing the range you want');
        rsRawPairs = [];
        if exist('Pairs2Use','var')
            Counts = feval(CCGFun,Res, uint32(Clu), BinSize, uint32(HalfBins),uint32(bUsePairs));
        else
            Counts = feval(CCGFun,Res, uint32(Clu), BinSize, uint32(HalfBins));
        end
    end
else
    if exist('Pairs2Use','var')
        Counts = feval(CCGFun,Res, uint32(Clu), BinSize, uint32(HalfBins),uint32(bUsePairs));
    else
        Counts = feval(CCGFun,Res, uint32(Clu), BinSize, uint32(HalfBins));
    end
end

% shape the results
nBins = 1+2*HalfBins;
% if there are no spikes in the top cluster, CCGEngine will produce a output the wrong size

Counts = double(reshape(Counts,[nBins nPresent nPresent]));
if nPresent<nGroups
    % extent array size with zeros
    Counts(nBins, nGroups, nGroups) = 0;
end
    
if nargout>=3
    Pairs = Included(ind(double(rsRawPairs) + 1));
end

% OK so we now have the bin counts.  Now we need to rescale it.

% remove bias due to edge effects - this should be vectorized
if isempty(Epochs)
    Bias = ones(nBins,1);
else
    nTerm = [HalfBins:-1:1 , 0.25 , 1:HalfBins];
	Bias = zeros(nBins,1);
    TotLen = 0;
	for e=1:size(Epochs,1)
        EpochLen = Epochs(e,2)-Epochs(e,1);
        EpochBias = clip(EpochLen - nTerm*BinSize,0,inf)*BinSize;
        Bias = Bias+EpochBias';
        TotLen = TotLen + EpochLen;
	end
    Bias = Bias/TotLen/BinSize;
end

if isempty(Epochs)
      Trange = max(Res) - min(Res); % total time
else
       Trange = sum(diff(Epochs,[],2));
end
t = 1000*(-HalfBins:HalfBins)*BinSize/SampleRate;

% count the number of spikes in each group:
for g=1:nGroups
	nSpikesPerGroup(g) = sum(Clu==g);
end;

% normalize each group
for g1=1:nGroups, for g2=g1:nGroups
	switch Normalization
		case 'hz'
			Factor = SampleRate / (BinSize * nSpikesPerGroup(g1));
			AxisUnit = '(Hz)';
		case 'hz2'
			Factor = SampleRate * SampleRate / (Trange*BinSize);	
			AxisUnit = '(Hz^2)';
		case 'count';
			Factor = 1;
			AxisUnit = '(Spikes)';
		case 'scale'
			Factor = Trange / (BinSize * nSpikesPerGroup(g1) * nSpikesPerGroup(g2));
			AxisUnit = '(Scaled)';
		otherwise
			warning(['Unknown Normalization method ', Normalization]);
	end;
% 	ccg(:,g1,g2) = flipud(Counts(:,g1,g2)) * Factor ./repmat(Bias,[1 nGroups,nGroups]); 
 	ccg(:,g1,g2) = (Counts(:,g1,g2)) * Factor ./Bias; 
	ccg(:,g2,g1) = flipud(Counts(:,g1,g2)) * Factor ./Bias; 
	% now plot, if there is no output argument
	if (nargout==0)
		FigureIndex = g1 + nGroups*(nGroups-g2);
		subplot(nGroups,nGroups,FigureIndex);		
	
		% plot graph
%		bar(1000*(-HalfBins:HalfBins)*BinSize/SampleRate, ccg(:,g1,g2));
		bar(t, ccg(:,g1,g2));

		% label y axis
		if g1==g2
%			ylabel(['ACG ', AxisUnit])	
    		FiringRate = SampleRate * nSpikesPerGroup(g1) / Trange;
%     		Ttitle = sprintf('%d (~%5.2fHz)',GSubset(g1),FiringRate);
% 			title(Ttitle);
    		xlabel('ms');
        else 
        %    set(gca, 'xtick', []);
        end
        if g1==1
            ylabel(sprintf('%d', GSubset(g2)));
%			ylabel(['CCG ', AxisUnit])
% 			Ttitle = sprintf('%d vs %d', GSubset(g1), GSubset(g2));
% 			title(Ttitle);
		end
        if g2==nGroups
    		Ttitle = sprintf('%d (~%5.2fHz)',GSubset(g1),FiringRate);
			title(Ttitle);
        end

		axis tight
	end
end,end;

% only make an output argument if its needed (to prevent command-line spew)
if (nargout>0)
	out = ccg;
end


