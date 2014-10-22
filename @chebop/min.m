function sol = min(A, L, f)
%MIN   Linear-quadratic optimization.
%   MIN(A, L) for linear chebops A and L solves the quadratic
%   minimization problem
%
%      Minimize    1/2 <u, A u>
%      Subject to  L u = f
%                  and all functional constraints attached to L.
%
%   The solution is returned as a chebmatrix whose components are
%   the minimizer u* and all adjoint variables. There is one adjoint
%   variable (a function) corresponding to u* and a scalar for each
%   constraint on L.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.


%%
% Check linearity.
try
    A = linop(A);
    L = linop(L);
catch ME
    if ( ME.identifier == 'CHEBFUN:CHEBOP:min:nonlinear' )
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
if ( A.diffOrder > L.diffOrder )
    error('CHEBFUN:CHEBOP:min:diffOrder', ...
        'diffOrder of objective A cannot be larger than that of constraints L.')
end

%%
% TODO: Check for cumsums, and throw an error if there are any.

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

%%
% Now get to work.

dom = A.domain;     % Domain of the problem
n = L.diffOrder;    % The maximum diffOrder of the problem
Z = operatorBlock.zeros(dom);

% Assembling the block operator is just this easy.
KKT = linop([ (A+A.')/2, L.'; L, Z ]);

% The boundary constraints take more work.
% First, construct vectors of operators
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

% Next, get the complementarity matrices of A and L.
CMA = complementarityMatrix(A, n);
CML = complementarityMatrix(L, n);

% Now loop through and add all the boundary conditions.
for k = 1:n
    % First BC is taking rows for the top-left corner, second
    % from the bottom-right.
    KKT = addbc(KKT, [CMA(k,1:n)*eval_a,       CML(k,1:n)*eval_a], 0);
    KKT = addbc(KKT, [CMA(n+k,n+1:end)*eval_b, CML(n+k,n+1:end)*eval_b], 0);
end

x = chebfun(@(x) x, dom);
sol = KKT \ [0*x; f];

end
