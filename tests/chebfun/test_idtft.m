function pass = test_idtft(pref)
% Test for IDTFT()

if ( nargin < 1 )
    pref = chebfunpref();
end

a = 1/2;
n = 53;
x = zeros(2*n+1, 1);
x(n+1:end)=a.^(0:n);
f = chebfun.dtft(x);
y = idtft(f);
pass = norm(x-y, inf) < 1e-13;
