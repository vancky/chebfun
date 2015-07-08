function x = idtft(X)
%DTFT   Inverse Discrete time Fourier transform of a chebfun X. 
%   x = DTFT(X) returns the inverse discrete time Fourier transform of the 
%   periodic chehbfun X, the domain of X must be [-pi, pi]. The 
%   returned x is a vector of  doubles. The definition of IDTFT is 
%   what is standard in signal processing:
%     x[n] = 1/(2*pi)int_{pi}^{pi} X(w)*exp(1i*w*n) dw

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isempty(X) )
    x = [];
    return
end
dom = domain(X);
a = dom(1);
b = dom(end);

if ( ~isa(X, 'chebfun') || ~isPeriodicTech(X) )
    error( 'CHEBFUN:IDTFT', 'X must be a periodic chebfun.')
end

if ( norm([a, b] - [-pi, pi], inf) > 1e-13 )
    error( 'CHEBFUN:IDTFT', 'domain of X must be [-pi, pi].')
end

%%
% get the coefficients and flip them due to the definition of DTFT:
x = flipud(get(X, 'coeffs'));
