function bol = isempty( N ) 
%ISEMPTY    Empty check for chebop2. 
% 
%  ISEMPTY( N ) returns true if N is an empty chebop2 object. Otherwise, it
%  returns false. 

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Every chebop2 must have an PDO: 
if ( ~isfield(N, 'op') || isempty( N.op ) )
    bol = 1; 
else
    bol = 0; 
end

end