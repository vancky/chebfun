function D = diff(disc, m)
%DIFF    Differentiation operator for COLLOCFOUR discretization.
%   D = DIFF(DISC) gives the matrix such that if v=D*u, then v=u', where u
%   is a COLLOCFOUR representation of a trigonometric polynomial.
%
%   DIFF(DISC, M) for positive integer M returns D^M (through a better
%   algorithm than multiplication).

%  Copyright 2014 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org/ for Chebfun information.

% Store information about domain and dimensions.
d = disc.domain;
n = disc.dimension;

if ( m == 0 )
    % Trivial case.
    D = eye(sum(n));
else
    numIntervals = disc.numIntervals;
    
    % Find the diagonal blocks.
    blocks = cell(numIntervals);
    for k = 1:numIntervals
        len = d(k+1) - d(k);
        % Scaled DIFFMATs.
        blocks{k} = collocFour.diffmat(n(k), m) * (2*pi/len)^m; 
    end
    
    % Assemble.
    D = blkdiag(blocks{:});
end

end