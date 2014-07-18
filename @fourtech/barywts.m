function w = barywts(n)
%BARYWTS   Barycentric weights for equally spaced points.
%   BARYWTS(N) returns the N barycentric weights for trigonometric polynomial 
%   interpolation on an equally spaced grid. The weights are normalised so 
%   that they have infinity norm equal to 1.
%
% See also BARY, FOURPTS.   

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% See ??
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ( n == 0 )                      % Special case (no points)
    w = [];
elseif ( n == 1 )                  % Special case (single point)
    w = 1;
else                               % General case
    w = ones(n, 1);
    w(2:2:end) = -1;
    w(1) = 1/2*w(1);
    w(end) = 1/2*w(end);
end

end
