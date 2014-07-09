function D = double(A)
%DOUBLE  Convert chebmatrix to double if possible.
%   If all the blocks in the chebmatrix are numeric, then DOUBLE(A) returns
%   the equivalent matrix with the same entries.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

try
    D = cell2mat(A.blocks);
catch
    error('Conversion to double not defined for this chebmatrix.')
end

end