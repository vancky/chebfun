function fx = trigBary(x, fvals, xk, vk)
%TRIGBARY   Trigonometric barycentric interpolation formula.
%   TRIGBARY(X, FVALS, XK, VK) uses the 2nd form barycentric formula with 
%   weights VK to evaluate an interpolant of the data {XK, FVALS(:,k)} at 
%   the points X. Note that XK and VK should be column vectors, and FVALS, 
%   XK, and VK should have the same length.
%
%   TRIGBARY(X, FVALS) assumes XK are equally spaced points in [-pi, pi) 
%   and VK are the corresponding barycentric weights.
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

% Default to equispaced nodes and barycentric weights:
if ( nargin < 3 )
    xk = linspace(-pi, pi, n+1);
    xk = xk(1:end-1); % discard the point x = pi;
end

if ( nargin < 4 )
    vk = fourtech.barywts(n);
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
