function varargout = chebfir(n, freqs, f)
%CHEBFIR   Chebfun filter desigin for FIR filters.
%
% See also CF.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

g = chebfun(@(x) f(1/pi*acos(x)), 'splitting', 'on');
dom = sort(cos(pi*freqs));
[p, err, status] = remez(g, n, 'domain', dom);
status.xk = sort(1/pi*acos(status.xk));
a = chebcoeffs(p);
c = [1/2*a(end:-1:2); a(1); 1/2*a(2:end);];
H = chebfun(c, 'periodic', 'coeffs');
varargout = {H, err, status};
end