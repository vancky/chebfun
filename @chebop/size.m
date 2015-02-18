function mn = size(N)
%SIZE   Size of a chebop.
%   SIZE(A) returns a vector [M, N], where M is the number of input
%   functions to the operator of the chebop A and N is the number of
%   its output functions.

% Evaluate N.op in order to determine the output size.
zeroFun = chebfun(0, N.domain);
u = repmat({zeroFun}, nargin(N.op), 1);
out = chebmatrix( N.op(u{:}) );
m = size(out, 1);

% Find out how many unknown variables N acts on.
n = numVars(N);

mn = [m n];

end
