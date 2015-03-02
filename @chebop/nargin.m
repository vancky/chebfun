function nIn = nargin(N)
%NARGIN   The number of input arguments to a CHEBOP .OP field.
%   NARGIN(N) returns the number of input arguments to N.OP, or zero if N.OP is
%   empty.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isempty(N.op) )
    nIn = 0;
elseif ( nargin(N.op) >= 0 )
    nIn = nargin(N.op);
else
    % Number of variables N acts on plus one for the independent variable.
    nIn = length(N.inDims) + 1;
end

end
