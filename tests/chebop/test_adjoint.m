function pass = test_adjoint(pref)
% Test file for the ADJOINT method.

if ( nargin == 0 )
    pref = cheboppref();
end
tol = 1e-10;

%% First-order:

% Construct an operator N = a1(x)u' + a0(x)u:
dom = [0 2*pi];
a1 = chebfun(@(x) exp(cos(x)), dom, 'trig');
a0 = chebfun(@(x) sin(10*x), dom, 'trig');
N = chebop(dom);
N.op = @(x,u) a1.*diff(u) + a0.*u;
N.bc = 'periodic';

% Compute the adjoint and its coefficients:
[Nstar, adjcoeffs] = adjoint(N);

% The exact adjoint is -a1(x)u' + (a0(x)-a1'(x))u. Let's compare the
% coefficients:
pass(1) = norm(adjcoeffs{1} - (-a1), inf) < tol;
pass(2) = norm(adjcoeffs{2} - (a0-diff(a1)), inf) < tol;

% Let's contruct the exact adjoint, and compare the evaluation against a 
% function:
f = chebfun(@(x) cos(20*x) + exp(sin(cos(x))), dom, 'trig');
Nstarex = chebop(@(x,u) -a1.*diff(u) + (a0-diff(a1)).*u, dom);
pass(3) = norm(Nstarex*f - Nstar*f, inf) < tol;


%% Second-order:

% Construct an operator N = a_2(x)*u'' + a_1(x)u' + a_0(x)u:
dom = [0 2*pi];
a2 = chebfun(@(x) cos(4*x), dom, 'trig');
a1 = chebfun(@(x) cos(cos(2*x)), dom, 'trig');
a0 = chebfun(@(x) exp(sin(x)), dom, 'trig');
N = chebop(dom);
N.op = @(x,u) a2.*diff(u,2) + a1.*diff(u) + a0.*u;
N.bc = 'periodic';

% Compute the adjoint:
[Nstar, adjcoeffs] = adjoint(N);

% The exact adjoint is N* = (a_2(x)u)'' - (a_1(x)u)' + a_0(x)u
% i.e., N* = a_2(x)u'' + (2a_2'(x)-a1(x))u' + (a_0(x)-a_1'(x)+a_2''(x))u.
% Let's compare the coefficients:
pass(4) = norm(adjcoeffs{1} - a2) < tol;
pass(5) = norm(adjcoeffs{2} - (2*diff(a2) - a1)) < tol;
pass(6) = norm(adjcoeffs{3} - (a0 - diff(a1) + diff(a2, 2))) < tol;

% Let's contruct the exact adjoint, and compare the evaluation against a 
% function:
f = chebfun(@(x) sin(cos(10*x)) + exp(sin(cos(x))), dom, 'trig');
Nstarex = chebop(@(x,u) diff(a2.*u, 2) - diff(a1.*u) + a0.*u, dom);
pass(7) = norm(Nstarex*f - Nstar*f, inf) < tol;

end
