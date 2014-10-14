function pass = test_transpose(pref)

if ( nargin == 0 )
    pref = cheboppref();
end

dom = [-2, -1];
x = chebfun(@(x) x, dom);
[Z, I, D, C, M] = linop.primitiveOperators(dom);
[z, e, s, r] = linop.primitiveFunctionals(dom);

% Make sure the transpose has no boundary conditions.
L = linop(D^2);
L = addConstraint(L, s, 2);
A = L.';
numConstraints = prod(size(A.constraint.functional));
pass(1) = (numConstraints == 0);

% Make sure the coefficients work properly.
L = linop(M(sin(x))*D^2);
A = L.';
Acoeff = toCoeff(A.blocks{1});
exact = chebmatrix([sin(x) 2*cos(x) -sin(x)]);
tol = 1e-12;
err = norm(Acoeff - exact);
pass(2) = err < tol;

% Make sure the transpose of the transpose is the original operator
% (sans boundary conditions).
AA = A.';
AAcoeff = toCoeff(AA.blocks{1});
Lcoeff = toCoeff(L.blocks{1});
err = norm(Lcoeff - AAcoeff)
pass(3) = (err == 0);

% Another check to make sure the transpose of the transpose is
% the same operator.
f = chebfun(@(x) sin(10*x) - exp(x), dom);
err = norm(L*f - AA*f)
pass(4) = (err == 0);

% Make sure an error is thrown if the linop has an integration operator.
L = linop(D + C);
pass(5) = 0;
try
    A = L.';
catch
    pass(5) = true;
end


end
