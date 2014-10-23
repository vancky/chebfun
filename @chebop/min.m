function sol = min(A, L, f, B, c)
%MIN   Linear-quadratic optimization.
%   MIN(A,L,f,B,c) for linear chebops A and L solves the quadratic
%   minimization problem
%
%      Minimize    1/2 <u, A u>
%      Subject to  L u = f
%             and  sum(B(k)*u) = c(k)
%                      for each chebfun B(k) and scalar c(k).
%
%   The solution is returned as a chebmatrix whose components are
%   the minimizer u* and all adjoint variables. There is one adjoint
%   variable (a function) corresponding to u* and a scalar for each
%   functional constraint.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.


%%
% Check linearity.
try
    A = linop(A);
    L = linop(L);
catch ME
    if ( strcmpi(ME.identifier, 'CHEBFUN:CHEBOP:min:nonlinear') )
        error('CHEBFUN:CHEBOP:min:nonlinear', ...
            'Nonlinear minimization is not supported.')
    else
        rethrow(ME)
    end
end

%%
% Test if either linop's size is any bigger than 1x1. If so, throw an error.
if ( (max(size(A)) > 1) || (max(size(L)) > 1) )
    error('CHEBFUN:CHEBOP:min:size', ...
        'Minimization for systems of equations is not supported.');
end

%%
% Check diffOrders.
% TOOD: Really A.diffOrder should be strictly less than that of L....
if ( A.diffOrder > L.diffOrder )
    error('CHEBFUN:CHEBOP:min:diffOrder', ...
        'diffOrder of objective A cannot be larger than that of constraints L.')
end

%%
% TODO: Check for cumsums, and throw an error if there are any.
% TODO: Check for complex numbers. Moreover, test to see if this works for
%       complex numbers.

%%
% Be sure the domains of the operators are identical.
if ( any(A.domain ~= L.domain) )
    error('CHEBFUN:CHEBOP:min:domains', ...
        'Domains of the two input operators must match.')
end

%%
% Check for constraints on A.
if ( length(A.constraint) > 0 )
    warning('CHEBFUN:CHEBOP:min:Aconstraint', ...
        'Ignoring constraints on objective operator A.')
end

%%
% Check for constraints on L. (At the moment these are not supported.)
if ( length(L.constraint) > 0 )
    warning('CHEBFUN:CHEBOP:min:Lconstraint', ...
        'Ignoring constraints on objective operator L (not yet supported).')
end

%%
% TODO: What to do about the RHS that might be specified in A or L?

dom = A.domain;     % Domain of the problem
x = chebfun(@(x) x, dom);

%%
% Set the RHS = 0 if it was not specified.
if ( nargin < 3 )
    % f not specified.
    f = 0*x;
end
if ( nargin < 4 )
    % B not specified.
    B = [];
end

%%
% TODO: Ensure B is the right size.

if ( nargin < 5 )
    % c not specified.
    c = zeros(size(B));
end

%%
% Assembling the block operator is just this easy.
% The functional constraints still need to be added (below).
[Z, I, D, C, M] = linop.primitiveOperators(dom);
[z, e, s, r]    = linop.primitiveFunctionals(dom);
KKT = [ (A+A.')/2, L.'; L, Z ];


%% BOUNDARY CONDITIONS FOR THE FUNCTIONAL CONSTRAINTS.
n = L.diffOrder;        % The maximum diffOrder of the problem
nbc = length(B);        % Number of functional constraints.
BCC = zeros(2*n, nbc);  % Matrix encoding coeffs of scalar Lagrange multipliers.
aug_funs   = [];
aug_fzero  = [];
aug_scalar = zeros(nbc);
aug_duals  = [];
aug_dzero  = [];
for k = 1:nbc
    func = B{k};
    BCC(  1:n,  k) = getDeltaCoeffsAt(func, dom(1), n);
    BCC(n+1:2*n,k) = getDeltaCoeffsAt(func, dom(end), n);

    % TODO: Treat deltas in the middle of domain, converting
    % to continuity conditions.
    funcSansDeltas = removeDeltas(func);
    aug_funs  = [aug_funs, funcSansDeltas];
    aug_fzero = [aug_fzero, 0*x];
    aug_duals = [aug_duals; r(func)]; % TODO: Must convert deltas to e() ops...
    aug_dzero = [aug_dzero; z];
end

% Add the functional constraints to the KKT system.
KKT = [ [KKT, [aug_funs; aug_fzero]];
         aug_duals, aug_dzero, aug_scalar ];


%% BOUNDARY CONDITIONS FOR THE ODE CONSTRAINT.
% Construct vectors of operators
%   [e, e*D, e*D*D, e*D*D*D, ...]
% that we can hit the coefficients against to get
% the proper boundary conditions. This is pretty slick.
D = operatorBlock.diff(dom);
eA = functionalBlock.feval(dom(1), dom);
eB = functionalBlock.feval(dom(end), dom);
eval_a = chebmatrix(eA);
eval_b = chebmatrix(eB);
for k = 1:n-1
    eval_a = [eval_a; eA*(D^k)];
    eval_b = [eval_b; eB*(D^k)];
end

% Get the complementarity matrices of A and L.
CMA = complementarityMatrix(A, n);
CML = complementarityMatrix(L, n);


%% NOW LOOP THROUGH AND ADD ALL THE BOUNDARY CONDITIONS.
KKT = linop(KKT);
for k = 1:n
    % Two BCs at a time (one for left endpoint, one for right).
    thisBC_a = [CMA(k,1:n)*eval_a,       CML(k,1:n)*eval_a];
    thisBC_b = [CMA(n+k,n+1:end)*eval_b, CML(n+k,n+1:end)*eval_b];
    for j = 1:nbc
        thisBC_a = [thisBC_a, BCC(k,  j)];
        thisBC_b = [thisBC_b, BCC(k+n,j)];
    end
    KKT = addbc(KKT, thisBC_a, 0);
    KKT = addbc(KKT, thisBC_b, 0);
end

% keyboard
sol = KKT \ [0*x; f; c];

end
