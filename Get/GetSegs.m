% [Segs Complete] = GetSegs(x, StartPoints, SegLen, IfNotComplete);
%
% extracts multiple segments from a time series x.  x may be a vector
% or a matrix, in which case time is on the first dimension, and channel
% is the second.
%
% StartPoints give the starting times of the segments.
% SegLen gives the length of the segments.  All segments must be
% the same length, and a rectangular array is returned
% (dim 1:time within segment, dim 2: segment number, dim 3:channel);
%
% IfNotComplete specifies what to do if the start or endpoints are outside
% of the range of x.  If a value is specified, any out-of-range points
% are given that number.  If [] is specified, incomplete segments are not
% returned.  A list of complete segments extracted is returned in Complete.
% Default value for IfNotComplete: NaN.

function [Segs, Complete] = GetSegs(x, StartPoints, SegLen, IfNotComplete);

if nargin<4
	IfNotComplete = NaN;
end

% make inputs column vectors.
if min(size(x)) == 1
	x = x(:);
end
StartPoints = StartPoints(:);

nSegs = length(StartPoints);
[xLen nDims] = size(x);

% all segments have the same length

% make index array for 1d array
IndMat1 = repmat(StartPoints(:)',[SegLen,1]) ...
			+ repmat((0:SegLen-1)', [1, nSegs]);
% find out of range points
OutOfRange1 = find(IndMat1<1 | IndMat1>xLen);
Complete = find(~any(IndMat1<1 | IndMat1>xLen));
IndMat1(OutOfRange1) = 1;

% now extend for multiple dimensions
IndMat = repmat(IndMat1, [1 1 nDims]) ...
		+ repmat(shiftdim((0:nDims-1)*xLen,-1), [SegLen, nSegs, 1]);

Segs = x(IndMat);

% replace out of range values
if size(IfNotComplete) == [1 1];
	for d=1:nDims
		Segs(OutOfRange1 + (d-1)*SegLen*nSegs) = IfNotComplete;
	end
elseif isempty(IfNotComplete)
	Segs(:,setdiff(1:nSegs, Complete),:) = [];
else
	error('IfNotComplete should be a scalar or []');
end

