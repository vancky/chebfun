function X = dtft(x)
%DTFT   Discrete time Fourier transform of a discrete signal x. 
%   X = DTFT(x) returns the discrete time Fourier transform of the signal
%   x[n] as a periodic chebfun on [-pi, pi]. x[n] must be a vector of 
%   doubles. The definition of DTFT used is the one standard in signal 
%   porocessing:
%     X(w) = sum_{n=-inf}^{n=inf} ( x[n]*exp(-1i*w*n) )
%   Notice the minus sign in the epxonential above.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isempty(x) )
    X = chebfun;
    return
end

if ( ~isvector(x) || ~isa(x, 'double') )
    error( 'CHEBFUN:DTFT', 'x must be a vector of doubles')
end

%%
% Make sure x is a column vector, and flip it due to the definition of
% DTFT:
x = flipud(x(:));
X = chebfun(x, [-pi, pi], 'coeffs', 'trig');
