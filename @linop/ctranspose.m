function B = ctranspose(A)
%CTRANSPOSE   Adjoint of a linear differential operator.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Test if the linop's size is any bigger than 1x1. If so, throw an error.
% We don't yet know exactly what to do for larger systems.
if ( max(size(A)) > 1 )
    error('CHEBFUN:LINOP:transpose:size', ...
        'TRANSPOSE only defined for 1x1 linops.');
end

dom = A.domain;     % The domain.
op = A.blocks{1};   % The operator block of interest.

% If the linop has any integration (CUMSUM) operators in it, then TOCOEFF
% gives an error. We catch that error and throw a more helpful one.
try
    coeffs = toCoeff(op); % Variable coefficients as a cell-array of chebfuns
catch ME
    if ( ME.identifier == 'CHEBFUN:BLOCKCOEFF:cumsum:notSupported' )
        % If the error was the CUMSUM error, then throw a more helpful one.
        error('CHEBFUN:LINOP:transpose:notSupported', ...
            'Adjoints of integration operators are not supported.');
    else
        % Otherwise, rethrow the original.
        rethrow(ME);
    end
end

[Z, I, D, C, M] = linop.primitiveOperators(dom);
m = op.diffOrder;       % Differential order of the operator
ad = Z;                 % Now create the formal adjoint
for k = 0:m             % Sum up terms
    ad = ad + (-1)^k * D^k * M(coeffs{m-k+1});
end

%%
% TODO: Need to extract the boundary coefficients before this works!
error('CHEBFUN:LINOP:ctranpose:notSupported', ...
    'Adjoints are not yet supported.')
U = BOUNDARY_COEFFICIENTS;

nu = rank(U);           % Number of conditions on u (op)
nw = 2*m - nu;          % Number of conditions on w (adjoint)

%%
% Extract the coefficients from the operator and evaluate them and their first
% k-1 derivatives at each of the endpoints. We will need these values later
% for the construction of the complementarity matrix.
zeta = [];
xi   = [];
for k = 1:m
    a_k = coeffs{m-k+1};           % The coeffs are indexed from highest order down.
    for p = 0:k-1                  % Now evaluate a_k and its diffs at endpoints.
        zeta(k,p+1) = a_k(dom(1)); % zeta(k,p+1) = (d^p a_k / dx^p)(a)
        xi(k,p+1) = a_k(dom(2));   %   xi(k,p+1) = (d^p a_k / dx^p)(b)
        a_k = diff(a_k);
    end
end

%%
% Next, build the complementarity matrix. The mathematics comes from
% integration by parts. See the 1960 paper by Greub and Rheinboldt for
% details.

% TODO: Vectorize this if possible.
MM = zeros(2*m);  % Complementarity matrix (to-be)
Z = zeros(m);
for k = 1:m
    A = zeros(m);
    B = zeros(m);
    for i = 0:k-1
        for j = 0:k-i-1
            if k-i-1 == 0,
                % Matlab doesn't like nchoosek(0,n).
                choose = 1;
            else
                choose = nchoosek(k-i-1,j);
            end
            A(i+1,j+1) = (-1)^(k-i)   * choose * zeta(k,k-i-j);
            B(i+1,j+1) = (-1)^(k-i-1) * choose * xi(k,k-i-j);
        end
    end
    MM = MM + [A, Z; Z, B];
end

%%
% Now the matrix W specifying the subspace that the adjoint BCs occupy
% comes from the product of a basis of NULL(U) and the complementarity
% matrix.
Q = null(U);
W = Q' * MM;
W = orth(W')';


end
