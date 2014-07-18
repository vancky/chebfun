function w = trigBaryWeights(x)
%BARYWEIGHTS   Barycentric weights for trigonometric interpolation.
%   W = TRIGBARYWEIGHTS(X) returns scaled barycentric weights for the points 
%   in the columm vector X. The weights are scaled such that norm(W, inf) == 1.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% [TODO]: Should this live in the trunk?

% Check inout dimension:
[n, m] = size(x);
if ( m > 1 )
    if ( n > 1 )
        error('CHEBFUN:trigBaryWeights:matrix', 'Input must be a vector.')
    else
        % Allow a row vector:
        n = m;
        x = x.';
    end
end

% [TODO]: How do we find the capacity for trig case?
% Capacity: 
% if ( isreal(x) )
%     C = 4/(max(x) - min(x));   % Capacity of interval.
% else
%     C = 1; % Scaling by capacity doesn't apply for complex nodes.
% end

% Compute the weights:
if ( (n < 2001) )              % For small n using matrices is faster.
   V = sin(1/2*bsxfun(@minus, x, x.'));
   V(1:n+1:end) = 1;
   VV = exp(sum(log(abs(V))));
   w = 1./(prod(sign(V)).*VV).';
   
else                           % For large n use a loop.
   w = ones(n,1);
   for j = 1:n
       v = sin(1/2*(x(j) - x)); v(j) = 1;
       vv = exp(sum(log(abs(v))));
       w(j) = 1./(prod(sign(v))*vv);
   end
end

% Scaling:
w = w./norm(w, inf);

end