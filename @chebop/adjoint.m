function [N, adjcoeffs] = adjoint(N)
%ADJOINT   Compute the adjoint of a linear CHEBOP.
%   N = ADJOINT(N), where N is a CHEBOP, returns the adjoint CHEBOP of N.
%
%   [N, adjcoeffs] = ADJOINT(N) also returns a CHEBMATRIX ADJCOEFFS which stores 
%   the (variables) coefficients of the adjoint. The indexation is as follows:
%      N = adjcoeffs{1}*u^(n) + adjcoeffs{2}*u^(n-1) + ... + adjcoeffs{n+1}*u
%
% See also LINOP/ADJOINT.

% Copyright 2015 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Tell CHEBOP/LINEARIZE to stop if it detects nonlinearity:
linCheck = true; 

% Linearize, thereby obtaining linearity information, a LINOP, and an input of
% the correct dimensions to pass to N:
[L, ~, isLinear, ~] = linearize(N, N.init, [], linCheck);

% We need the entire operator (including BCs) to be linear:
assert(all(isLinear), 'CHEBFUN:CHEBOP:adjoint:nonlinear', ...
    ['The input operator appears to be nonlinear.\n', ...
    'ADJOINT supports only linear CHEBOP instances.']);

% Get the value of the highest derivative:
n = L.diffOrder;

% ADJOINT is supported only for periodic boundary conditions for the moment.
% [TODO]: Support non-periodic boundary conditions.
if ( ( ~isa(N.bc, 'char') || ~strcmpi(N.bc,'periodic') ) && ( n > 0 ) ) 
    error('CHEBFUN:CHEBOP:adjoint:nonperiodic', ...
        'ADJOINT only supports periodic boundary conditions for the moment.');
end

% Extract the blocks of the LINOP:
blocks = L.blocks;

% ADJOINT doesn't support system of equations for the moment.
% [TODO]: Support system of equations.
if ( max(size(blocks,2)) > 1 )
    error('CHEBFUN:LINOP:adjoint:system', ...
        'ADJOINT doesn''t support system of equations for the moment.');
end

% Get the first block and the domain:
block = blocks{1};
dom = L.domain;

% Get the coefficients:
if ( n == 0 )
    % This is a multiplication operator, nothing to do here:
    if ( nargout > 1 )
        adjcoeffs = conj(block);
    end
    return
else
    coeffs = toCoeff(block);
end

% Compute the coefficients of the adjoint:
adjcoeffs = cell(n+1,1);
for k = 0:n
    adjcoeffs{k+1} = chebfun('0', dom);
end
for k = 0:n
    for l = 0:k
        adjcoeffs{n+1-l} = adjcoeffs{n+1-l} + ...
            (-1)^k*nchoosek(k,l)*conj(diff(coeffs{n+1-k}, k-l));
    end
end

% Construct a CHEBOP from these new coefficients:
N = chebop(@(x, u) coeffs2func(u, adjcoeffs), dom);
N.bc = 'periodic';

if ( nargout > 1 )
    % Store the coefficients in a CHEBMATRIX:
    adjcoeffs = chebmatrix(adjcoeffs);
end

end

function out = coeffs2func(u, coeffs)
%COEFFS2FUNC   Output a function representing the operator field of the adjoint.

out = 0;
n = length(coeffs);
for k = 1:n
    out = out + coeffs{k}.*diff(u, n-k);
end

end
