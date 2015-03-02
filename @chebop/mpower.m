function C = mpower(A, pow)
%^   Repeated application of a chebop.
%   A^P, where A is a CHEBOP and P a nonnegative integer, returns the CHEBOP
%   representing P-fold application of A. The chebop A must be scalar; i.e.
%   it must act on only one function.
%
% See also CHEBOP/MTIME.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Check second argument is valid.
if ( pow ~= round(pow) || pow < 0 )
    error('CHEBFUN:CHEBOP:mpower:badPower', ...
        'Power must be a positive integer.')
end

% Construct a CHEBOP for repeated application.
C = A;
for k = 2:pow
    % Any dimension errors will be caught in MTIMES.
    C = C*A;
end

end
