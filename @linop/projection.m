function v = projection(A,u)
%PROJECTION  Project function into the kernel.
% V = PROJECTION(A,U) returns a chebfun V such that A*V=0, for linop A. The
% difference V-U is as of low degree as possible--that is, size(A,1)-1. 

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% The minimal norm solution is not so easy, so we choose the minimum degree
% Chebyshev polynomial.

m = size(A,1);
T = chebpoly(0:m-1,domain(u));   % Chebyshev polynomials
c = double(A*T) \ double(A*u);   % solve for v in the T basis
v = u - T*c; 

end