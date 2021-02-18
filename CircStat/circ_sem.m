function s = circ_sem(x)

%circ_sem - Compute Circular standard error of the mean (SEM).
%
%  USAGE
%
%    s = circ_(x)
%
%    s              vector or matrix over which the SEM should be computed

% Pascale July 2017

if nargin < 1,
  error('Incorrect number of parameters (type ''help <a href="matlab:help sem">sem</a>'' for details).');
end

if ~isdmatrix(x) && ~isdvector(x),
  error('Incorrect input - use vector or matrix (type ''help <a href="matlab:help sem">sem</a>'' for details).');
end

if any(size(x)==1), x = x(:); end

s = circ_std(x)/sqrt(size(x,1));