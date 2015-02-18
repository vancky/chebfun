function out = vertcat(varargin)
%VERTCAT   Vertically concatenate CHEBOP objects.
%   Z = [A; B; C; ...] vertically combines the chebmatrices, chebops,
%   chebfuns, and scalars given in the call, if their column sizes are
%   compatible.

out = chebmatrix(varargin.');

end
