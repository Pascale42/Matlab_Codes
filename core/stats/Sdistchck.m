function [errorcode,varargout] = Sdistchck(nparms,varargin)
%DISTCHCK Checks the argument list for the probability functions.

%   B.A. Jones  1-22-93

errorcode = 0;
n = nargout-1;
varargout = cell(1,n);

if nparms == 1
    varargout{1} = varargin{1};
    return;
end

% Get size of each input, check for scalars, copy to output
sz = cell(1,n);
isscalar = logical(zeros(1,n));
for j=1:n
   s = size(varargin{j});
   sz{j} = s;
   isscalar(j) = (prod(s) == 1);
   varargout{j} = varargin{j};
end

% Done if all inputs are scalars.  Otherwise fetch their common size.
if (all(isscalar)), return; end
t = sz(~isscalar);
size1 = t{1};

% Scalars receive this size.  Other arrays must have the proper size.
for j=1:n
   sizej = sz{j};
   if (isscalar(j))
      t = zeros(size1);
      t(:) = varargin{j};
      varargout{j} = t;
   elseif (~isequal(sizej,size1))
      errorcode = 1;
      return;
   end
end
