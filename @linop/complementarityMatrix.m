function MM = complementarityMatrix(L, N)
%COMPLEMENTARITYMATRIX   Complementarity matrix for a linear operator.
%   COMPLEMENTARITYMATRIX(L,N) returns the 2N-by-2N matrix of
%   complementarity between the boundary conditions of a differential
%   operator L and its adjoint.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Test if the linop's size is any bigger than 1x1. If so, throw an error.
if ( max(size(L)) > 1 )
    error('CHEBFUN:LINOP:complementarityMatrix:size', ...
        'COMPLEMENTARITYMATRIX defined only for 1x1 linops.');
end

dom = L.domain;     % The domain.
op = L.blocks{1};   % The operator block of interest.
m = op.diffOrder;   % Differential order of the operator

if ( nargin < 2 )
    % Default to diffOrder size.
    N = m;
elseif ( N < m )
    % The size of the matrix must be at least m.
    error('CHEBFUN:LINOP:complementarityMatrix:dimension', ...
        'Size of complementarity matrix must be at least as large as diffOrder.')
end

% If the linop has any integration (CUMSUM) operators in it, then TOCOEFF
% gives an error. We catch that error and throw a more helpful one.
try
    coeffs = toCoeff(op); % Variable coefficients as a cell-array of chebfuns
catch ME
    if ( ME.identifier == 'CHEBFUN:BLOCKCOEFF:cumsum:notSupported' )
        % If the error was the CUMSUM error, then throw a more helpful one.
        error('CHEBFUN:LINOP:complementarityMatrix:notSupported', ...
            'Complementarity of integration operators is not supported.');
    else
        % Otherwise, rethrow the original.
        rethrow(ME);
    end
end

%%
% Extract the coefficients from the operator and evaluate them and their first
% k-1 derivatives at each of the endpoints. We will need these values for the
% construction of the complementarity matrix.
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
% integration by parts. See either Hrothgar's thesis-to-be or the the 1960
% paper by Greub and Rheinboldt for details.

% TODO: Vectorize this if possible.
MM = zeros(2*N);    % Complementarity matrix, 2N-by-2N
Z = zeros(N);       % Zeros matrix, N-by-N
for k = 1:m         % Compute entries up to the diffOrder
    A = zeros(N);   % But the submatrices might be larger than m
    B = zeros(N);
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


end
