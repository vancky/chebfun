function fx = trigBary(x, fvals, xk, dom)
%TRIGBARY   Trigonometric barycentric interpolation formula.
%   TRIGBARY(X, FVALS, XK, dom) uses the 2nd form barycentric formula 
%   to evaluate an interpolant of the data {XK, FVALS(:,k)}
%   at the points X. The interpolant is supposed to live on the domain 
%   specified in dom. Note that XK should be column vector, and 
%   FVALS should have the same length.
%
%   TRIGBARY(X, FVALS) assumes XK are equally spaced points in [-pi, pi).
%
%   If size(FVALS, 2) > 1 then TRIGBARY(X, FVALS) returns values in the form
%   [F_1(X), F_2(X), ...], where size(F_k(X)) = size(X) for each k.
%
% See also BARY

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Parse inputs:
[n, m] = size(fvals);
sizex = size(x);
ndimsx = ndims(x);

if ( (m > 1) && (ndimsx > 2) )
    error('CHEBFUN:trigbary:evalArrayAtNDArray', ...
        ['TRIGBARY does not support evaluation of vectors of polynomials at ' ...
         'inputs with more than two dimensions.']);
end

% If a domain is provided:
if ( nargin == 4 )
    % Map the points to [-pi, pi]
    a = dom(1);
    b = dom(2);
    xk = pi/(b-a)*(2*xk-a-b);
    x  = pi/(b-a)*(2*x-a-b);
    % Compute the weights
end

if ( any((xk > pi) | (xk < -pi)) )
     error('CHEBFUN:trigbary:invalidNodes', 'nodes XK must lie within the domain');
end

% Default to equispaced nodes in [-pi, pi) and barycentric weights:
if ( nargin < 3 )
    % This is more efficient than calling fourtech.fourpts.
    xk = linspace(-pi, pi, n+1);
    % Discard the point x = pi;
    xk = xk(1:end-1);
    % Default weights:
    vk = fourtech.barywts(n);
else
    % Compute the weights:
    vk = trigBarywts(xk);
end

if ( nargin < 4 )
    
end

if ( ~all(sizex) )
    fx = x;
    return
end

% Check that input is a column vector:
if ( (ndimsx > 2) || (sizex(2) > 1) )
    x = x(:);
end

% The function is a constant.
if ( n == 1 )
    fx = repmat(fvals, length(x), 1);
    return
end

% The function is NaN.
if ( any(isnan(fvals)) )
    fx = NaN(length(x), m);
    return
end

% Choose the appropriate function based on the length of the values to be
% interpolated:
if ( rem(n, 2) == 0 )
    s = cot(sum(xk)/2);
    ctsc = @(x) cot(x) + s;
else
    ctsc = @(x) csc(x);
end

% The main loop:
% Initialise:
num = zeros(size(x, 1), m);
denom = num;

% Loop:
for j = 1:length(xk),
    tmp = vk(j) * ctsc((x - xk(j))/2);
    num = num + bsxfun(@times, tmp, fvals(j,:));
    denom = bsxfun(@plus, denom, tmp);
end
fx = num ./ denom;

% Try to clean up NaNs:
for k = find(isnan(fx(:,1)))'       % (Transpose as Matlab loops over columns)
    index = find(x(k) == xk, 1);    % Find the corresponding node
    if ( ~isempty(index) )
        fx(k,:) = fvals(index,:);   % Set entry/entries to the correct value
    end
end

% Reshape the output if possible:
if ( (m == 1) && ( (ndimsx > 2) || (sizex(2) > 1) ) )
    fx = reshape(fx, sizex);
elseif ( (m > 1) && ( (ndimsx == 2) || (sizex(2) > 1) ) )
    fx = reshape(fx, sizex(1), m*numel(x)/sizex(1));
end

end