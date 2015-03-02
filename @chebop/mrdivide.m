function C = mrdivide(A, B)
%/    CHEBOP right divide.
%   C = A/B, where B is scalar, returns a CHEBOP C representing scalar
%   division of the original operator.
%
% See also CHEBOP/MTIMES.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isscalar(B) )
    C = (1/B)*A;

else
    error('CHEBFUN:CHEBOP:mrdivide:nonScalar', ...
        'CHEBOP / DOUBLE division is defined only for scalars.');

end
