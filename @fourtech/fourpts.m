function x = fourpts(n, dom)
%FOURPTS   Fourier points in [-1, 1).
%   FOURPTS(N) returns N equispaced points in [-1, 1).
%
% See also CHEBPTS, LEGPTS, JACPTS, LAGPTS, and HERMPTS.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

if ( nargin == 1 )
    dom = [-1,1];
end

x = linspace(dom(1), dom(end), n+1).';
x(end) = [];

end
