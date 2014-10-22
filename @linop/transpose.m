function B = transpose(A)
%TRANSPOSE   Formal adjoint of a linear differential operator.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Test if the linop's size is any bigger than 1x1. If so, throw an error.
% We don't yet know exactly what to do for larger systems.
if ( max(size(A)) > 1 )
    error('CHEBFUN:LINOP:transpose:size', ...
        'TRANSPOSE defined only for 1x1 linops.');
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
ad = Z;                 % Create the formal adjoint
for k = 0:m             % Sum up terms
    ad = ad + (-1)^k * D^k * M(coeffs{m-k+1});
end

% Construct a linop out of the adjoint operatorBlock.
B = linop(ad);

end
