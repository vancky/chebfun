function B = ctranspose(L)
%CTRANSPOSE   Adjoint of a linear differential operator.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Test if the linop's size is any bigger than 1x1. If so, throw an error.
% We don't yet know exactly what to do for larger systems.
if ( max(size(L)) > 1 )
    error('CHEBFUN:LINOP:ctranspose:size', ...
        'TRANSPOSE defined only for 1x1 linops.');
end

dom = L.domain;     % The domain.
op = L.blocks{1};   % The operator block of interest.

% If the linop has any integration (CUMSUM) operators in it, then TOCOEFF
% gives an error. We catch that error and throw a more helpful one.
try
    coeffs = toCoeff(op); % Variable coefficients as a cell-array of chebfuns
catch ME
    if ( ME.identifier == 'CHEBFUN:BLOCKCOEFF:cumsum:notSupported' )
        % If the error was the CUMSUM error, then throw a more helpful one.
        error('CHEBFUN:LINOP:ctranspose:notSupported', ...
            'Adjoints of integration operators are not supported.');
    else
        % Otherwise, rethrow the original.
        rethrow(ME);
    end
end

%%
% TODO: Need to extract the boundary coefficients before this works!
% TODO: Will we need to check to see if the BCs are homogeneous?
%%
% But for now we can return the adjoint of an operator with no BCs.
if ( length(L.constraint) == 0 )
    % It is possible to get them using the complementarity matrix by
    % specifying the following BC encoding matrix, but it is super inefficient
    % to do so because if L has no constraints we already know what the
    % constraints on the adjoint will be.
    %   U = zeros(0, 2*m);

    m = op.diffOrder; % Differential order of the operator
    B = L.';          % The formal adjoint

    D = operatorBlock.diff(dom);
    e = functionalBlock.eval(dom);

    for k = 0:m-1
        % Add each boundary condition.
        B = addbc(B, e(dom(1))*D^k, 0);
        B = addbc(B, e(dom(end))*D^k, 0);
    end

    return
else
    error('CHEBFUN:LINOP:ctranpose:notSupported', ...
        'Adjoints are not yet supported for operators with boundary conditions.')
end

%%
% Now the matrix W specifying the subspace that the adjoint BCs occupy
% comes from the product of a basis of NULL(U) and the complementarity
% matrix.
MM = complementarityMatrix(L);
Q = null(U);
W = Q' * MM;
W = orth(W')';  % Do we want to do this?


end
