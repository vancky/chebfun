function u = enforceBC(N,u)
%ENFORCEBC  Perturb function to satisfy boundary conditions. 
% V = ENFORCEBC(A,U) returns a chebfun V such that all the boundary
% conditions attached to N (lbc/rbc/bc properties) hold true for V. The
% perturbation is as of low polynomial degree as possible.
%
% This process usually works only if U is fairly close to satisfying them.
% If it fails, the original U is returned after a warning.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

u0 = u;
for newton = 1:5
    [L,val] = BCjacobian(N,u);       % linearize the BC at u
    m = size(L,1);
    T = chebpoly(0:m-1,domain(u));   % Chebyshev polynomials
    c = double(L*T) \ val;           % solve for delta in the T basis
    delta = -(T*c);
    u = u + delta;
    if (norm(delta) < 1e-13)
        return
    end
end

% Allow some tolerance here. 
if (norm(delta) > 1e-10)
    warning('BCs cannot be enforced correctly.')
    u = u0;
end

end

function [BCop,BCval] = BCjacobian(N,u)
BCop = [];
BCval = [];

uAD = adchebfun(u);
dom = domain(u);
x = chebfun(@(x) x,dom);

if ~isempty(N.bc)
    result = N.bc(x,uAD);
    for k = 1:length(result)
        r = result{k};
        BCop = [ BCop; linop(r.jacobian) ];
        BCval = [ BCval; double(r.func) ];
    end
end

[~,E] = linop.primitiveFunctionals(dom);

if ~isempty(N.lbc)
    e = E( dom(1) );
    result = N.lbc(uAD);
    for k = 1:length(result)
        r = result{k};
        BCop = [ BCop; e*linop(r.jacobian) ];
        BCval = [ BCval; r.func(dom(1)) ];
    end
end

if ~isempty(N.rbc)
    e = E( dom(2) );
    result = N.rbc(uAD);
    for k = 1:length(result)
        r = result{k};
        BCop = [ BCop; e*linop(r.jacobian) ];
        BCval = [ BCval; r.func(dom(2)) ];
    end
end

end
