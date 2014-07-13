function vscl = vscale(f) 
%VSCALE   Vertical scale of a CHEBFUN2.
% 
% VSCL = VSCALE(F) returns the vertial scale of a CHEBFUN2 as determined
% by evaluating on a coarse Chebyshev tensor-product grid. 

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% TODO: Should this also be taking the maximum along the edges when we are
% evaluating at 1st kind grids. 

% If f is an empty Chebfun2, VSCL = 0: 
if ( isempty(f) ) 
    vscl = 0; 
    return
end

% Get the degree of the CHEBFUN2:[m, n] = length(f); 
[m, n] = length( f ); 
%vscl = max( abs( f.pivotValues(1) ) ); 

% If F is of low degree, then oversample: 
m = max(m, 9); 
n = max(n, 9); 

% Calculate values on a tensor grid: 
X = rot90( chebcoeffs2( f ), 2 ); 
X(size(X,1)+1:m, :) = 0; X(:, size(X,2)+1) = 0;  % pad
X = rot90(X, 2); 
vals = chebfun2.coeffs2vals( X ); 

% Take the absolute maximum: 
vscl = max(abs(vals(:))); 

end