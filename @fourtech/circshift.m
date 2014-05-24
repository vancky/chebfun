function f = circshift(f,a)
%CIRCSHIFT   Circular (periodic) shift of a FOURTECH
%   G = CIRCSHIFT(F,A) circularly (periodically) shift F by A, where A is 
%   a real-valued number.  The result is a new FOURTECH G that is equal
%   to G = F(X-A). If F is an array-valued FOURTECH then each column of F
%   is shifted by A.
%
% See also FOURTECH.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org for Chebfun information.

[n,m] = size(f);

if numel(a) > 1
    error('CHEBFUN:FOURTECH:circshift', ...
          'Only scalar-valued shifts are allowed')
end
% Trivial case
if ( a == 0 )    
    return
end

% The FOURTECH object is defined as 
%   If N is odd
%       F(x) = C(1)*z^(N-1)/2 + C(2)*z^((N-1)/2-1) + ... + C(N)*z^(-(N-1)/2)
%   If N is even
%       F(x) = C(1)*z^(N/2-1) + C(2)*z^(N/2-2) + ... + C(N-1)*z^(-N/2-1) +
%                  1/2*C(N)*(z^(N/2) + z^(-N/2))
%   where z = exp(1i*pi*x) and -1 <= x <= 1. 

% So a shift by a only requires multiplying the kth coefficients by
% (exp(-1i*pi*a))^k

if mod(n,2) 
    even_odd_fix = (exp(-1i*pi*a)).^((n-1)/2:-1:-(n-1)/2);
else
    even_odd_fix = (exp(-1i*pi*a)).^((n/2-1):-1:(-n/2));
end
f.coeffs = bsxfun(@times,f.coeffs,even_odd_fix.');
f.values = f.coeffs2vals(f.coeffs);
if ~isreal(a)
    f.isReal = false(1,size(f.coeffs,2));
end

f.values(:,f.isReal) = real(f.values(:,f.isReal));

end