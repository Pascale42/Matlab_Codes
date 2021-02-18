function [SegsPow, f, t] = SegsPower(Segs,  mtFunction, varargin)

nSegs = size(Segs,3);
SegsPow =[];
for i=1:nSegs
    nArgs =nargin-1;
    fargs = cell(1,nArgs);
    fargs{1}= sq(Segs(:,:,i));
    for argi=2:nArgs
        fargs{argi} = varargin{argi-1};
    end
    try  
        [y, f, t] = feval(mtFunction, fargs);
    catch
        [y, f] = feval(mtFunction, fargs);
    end
    
    nDims = ndims(y);
    SegsPow = cat(nDims+1, SegsPow, y);
end