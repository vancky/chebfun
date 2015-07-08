function pass = test_dtft(pref)
% Test for DTFT()

if ( nargin < 1 )
    pref = chebfunpref();
end


a = 1/2;
n = 53;
x = zeros(2*n+1, 1);
x(n+1:end)=a.^(0:n);
f = chebfun.dtft(x);
F = chebfun(@(w) 1./(1-a*exp(-1i*w)), [-pi, pi], 'trig');
pass = norm(f-F, inf) < 1e-13;
