function BV = boundaryValues(F, N)
%BOUNDARYVALUES   Boundary values of a chefun.
%   BOUNDARYVALUES(F,N) computes the boundary values of a chebfun F and its
%   first N-1 derivatives, returning the result in a vector of the form
%       b = [F(a), F'(a), F''(a), ..., F(b), F'(b), F''(b), ...].'

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( N < 1 )
    BV = [];
    return
end

dom = F.domain;
a   = dom(1);
b   = dom(end);

aa = [];
bb = [];
aa(1) = feval(F, a);
bb(1) = feval(F, b);

% TODO: Make this (possibly) more efficient by checking the length of F and
% don't compute derivatives that are guaranteed to be zero.
for k = 2:N
    F = diff(F);
    aa(k) = feval(F, a);
    bb(k) = feval(F, b);
end

BV = [aa(:); bb(:)];

end
