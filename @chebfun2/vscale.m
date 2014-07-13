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

% Calculate absolute value of first GE pivot: 
maxCols = vscale( f.cols(:, 1) ); 
maxRows = vscale( f.rows(:, 1) );
pivValue = abs( f.pivotValues(1) ); 
vscl = maxCols * maxRows ./ pivValue; 

end