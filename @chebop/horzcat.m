function out = horzcat(varargin)
%HORZCAT   Horizontally concatenate CHEBOP objects.
%   Z = [A, B, C, ...] horizontally combines the chebmatrices, operator blocks,
%   chebfuns, and scalars given in the call, if their column sizes are
%   compatible.

out = chebmatrix(varargin);

end
