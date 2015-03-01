function varargout = firpm(n, freqs, f, varargin)
%FIRPM   FIR filter design in chebfun. 
%   H = FIRPM(N, FREQS, F) computes the best minimax degree N filter of length
%   2N+1 to match the desired response F specified by bands in FREQS.
%
%   N is the degree of the trigonometric polynomial used to design the filter.
%   The number of filter coefficients is always 2N+1.
%
%   FREQS is an even length vector of monotonically increasing frequencies 
%   in [0, 1] with FREQS(1) = 0 and FREQS(end) = 1. The frequency interval 
%   specified by [FREQS(k), FREQS(k+1)] for k odd is considered a don't care 
%   frequency band.
%
%   F must be a chebfun that lives on [0, 1] and specifies the ideal 
%   frequency response of the filter.
%
%   [...] = FIRPM(..., 'tol', TOL) uses the value TOL as the termination
%   tolerance on the increase of the levelled error.
%
%   [...] = FIRPM(..., 'display', 'iter') displays output at each iteration.
%
%   [...] = FIRPM(..., 'maxiter', MAXITER) sets the maximum number of allowable
%   iterations to MAXITER.
%
%   [...] = FIRPM(..., 'plotfcns', 'error') plots the error after each iteration
%   while the algorithm executes.
%
%   [P, ERR] = FIRPM(...) returns the maximum error ERR.
%
%   [P, ERR, STATUS] = FIRPM(...)  also return a structure array STATUS 
%   with the following fields:
%      STATUS.DELTA  - Obtained tolerance.
%      STATUS.ITER   - Number of iterations performed.
%      STATUS.DIFFX  - Maximum correction in last trial reference.
%      STATUS.XK     - Last trial reference on which the error equioscillates.
%
%
% References:
%
%   [1] Javed, M. and Trefethen, L. N.  "Remez and CF approximations of 
%    Periodic Functions". In preparation.
%
% See also REMEZ, TRIGCF, CF.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( nargin < 3 )
    error('CHEBFUN:CHEBFUN:firmp:nargin', ...
        'FIRPM should have at least three input arguments.');
end

if ( isempty(f) )
    varargout = {f};
    return;
end

if ( ~isreal(f) )
    error('CHEBFUN:CHEBFUN:firmpm:real', ...
        'FIRPM only supports real-valued functions.');
end

if ( numColumns(f) > 1 )
    error('CHEBFUN:CHEBFUN:firmpm:quasi', ...
        'FIRPM does not currently support quasimatrices.');
end

if ( issing(f) )
    error('CHEBFUN:CHEBFUN:firmpm:singularFunction', ...
        'FIRPM does not support functions with singularities.');
end

if ( isdelta(f) )
    error('CHEBFUN:CHEBFUN:firmpm:deltaFunctions', ...
        'Desired frequency resoonse must not have any delta functions.');
end

dom = f.domain([1, end]);
normf = norm(f);
if ( any((abs(dom - [0, 1]) > 100*eps )) )
    error('CHEBFUN:CHEBFUN:firmpm:Domain', ...
        'Desired frequency response must be defined on [0, 1].');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Wrapper to call the pmRemez Algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ( isa(f, 'double') )
    f = chebfun.interp1(freqs, f, 'linear', [0, 1]);
end

% Map the problem from the circle onto [-1, 1]:
g = chebfun(@(x) feval(f, 1/pi*acos(x)), 'splitting', 'on');


bands = sort(cos(pi*freqs));
[p, err, status] = pmRemez(g, n, bands);
status.xk = sort(1/pi*acos(status.xk));
a = chebcoeffs(p);
c = [1/2*a(end:-1:2); a(1); 1/2*a(2:end);];
H = chebfun(c, 'trig', 'coeffs');
varargout = {H, err, status};
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PM version of the Remez Algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = pmRemez(f, n, bands)
% Parse the inputs.
[n, N, opts] = parseInputs(f, n);
normf = norm(f);
% Initial values for some parameters.
iter = 0;       % Iteration count.
delta = normf;  % Value for stopping criterion.
deltamin = inf; % Minimum error encountered.
diffx = 1;      % Maximum correction to trial reference.

% Compute an initial reference set to start the algorithm.
xk = getInitialReference(f, n, bands);
xo = xk;

% Print header for text output display if requested.
if ( opts.displayIter )
    disp('It.     Max(|Error|)       |ErrorRef|      Delta ErrorRef      Delta Ref')
end

% Run the main algorithm.
while ( (delta/normf > opts.tol) && (iter < opts.maxIter) && (diffx > 0) )
    fk = feval(f, xk);     % Evaluate on the exchange set.
    w = baryWeights(xk);   % Barycentric weights for exchange set.

    [p, h] = computeTrialFunctionPolynomial(fk, xk, w, n, N, bands);
    
    % Perturb exactly-zero values of the levelled error.
    if ( h == 0 )
        h = 1e-19;
    end

    % Update the exchange set using the Remez algorithm with full exchange.
    [xk, err, err_handle] = exchange(xk, h, 2, f, p, N + 2, bands);

    % If overshoot, recompute with one-point exchange.
    if ( err/normf > 1e5 )
        [xk, err, err_handle] = exchange(xo, h, 1, f, p, N + 2, bands);
    end

    % Update max. correction to trial reference and stopping criterion.
    diffx = max(abs(xo - xk));
    delta = err - abs(h);

    % Store approximation with minimum norm.
    if ( delta < deltamin )
        pmin = p;        
        errmin = err;
        xkmin = xk;
        deltamin = delta;
    end

    % Display diagnostic information as requested.
    if ( opts.plotIter )
        doPlotIter(xo, xk, err_handle, dom);
    end

    if ( opts.displayIter )
        doDisplayIter(iter, err, h, delta, normf, diffx);
    end

    xo = xk;
    iter = iter + 1;
end

% Take best results of all the iterations we ran.
p = pmin;
err = errmin;
xk = xkmin;
delta = deltamin;

% Warn the user if we failed to converge.
if ( delta/normf > opts.tol )
    warning('CHEBFUN:CHEBFUN:pmRemez:convergence', ...
        ['PM-Remez algorithm did not converge after ', num2str(iter), ...
         ' iterations to the tolerance ', num2str(opts.tol), '.']);
end

% Form the outputs.
status.delta = delta/normf;
status.iter = iter;
status.diffx = diffx;
status.xk = xk;

varargout = {p, err, status};

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input parsing.

function [n, N, opts] = parseInputs(f, n)

varargin = {};
N = n;

% Parse name-value option pairs.
opts.tol = 1e-16*(N^2 + 10); % Relative tolerance for deciding convergence.
opts.maxIter = 40;           % Maximum number of allowable iterations.
opts.displayIter = false;    % Print output after each iteration.
opts.plotIter = false;       % Plot approximation at each iteration.

for k = 1:2:length(varargin)
    if ( strcmpi('tol', varargin{k}) )
        opts.tol = varargin{k+1};
    elseif ( strcmpi('maxiter', varargin{k}) )
        opts.maxIter = varargin{k+1};
    elseif ( strcmpi('display', varargin{k}) )
        opts.displayIter = true;
    elseif ( strcmpi('plotfcns', varargin{k}) )
        opts.plotIter = true;
    else
        error('CHEBFUN:CHEBFUN:pmRemez:badInput', ...
            'Unrecognized sequence of input parameters.')
    end
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions implementing the core part of the algorithm.

function xk = getInitialReference(f, N, bands)

% In the polynomial case just use the Chebyshev points:
xk = chebpts(N + 2, f.domain([1, end]));

if ( length(bands) == 4 )
    L = bands(2)-bands(1) + bands(4)-bands(3);
    n1 = (N+2)*(bands(2)-bands(1))/L;
    n2 = (N+2)*(bands(4)-bands(3))/L;
    np = ceil(n1);
    ns = floor(n2);
    if ( np + ns > N + 2 )
        np = np - 1;
    end        
    
    if ( np + ns < N + 2 )
        ns = ns + 1;
    end
    
    xk1 = chebpts( np, bands(1:2));
    xk2 = chebpts( ns, bands(3:4));
    xk = [xk1; xk2];
    if ( length(xk) ~= N + 2 )
        error( 'length is not N+2' )
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [p, h] = computeTrialFunctionPolynomial(fk, xk, w, m, N, bands)

% Vector of alternating signs.
sigma = ones(N + 2, 1);
sigma(2:2:end) = -1;

h = (w'*fk) / (w'*sigma);                          % Levelled reference error.
pk = (fk - h*sigma);                               % Vals. of r*q in reference.

% Trial polynomial.
p = chebfun(@(x) bary(x, pk, xk, w), [-1, 1], m + 1);

end

function [xk, norme, err_handle, flag] = exchange(xk, h, method, f, p, Npts, bands)
%EXCHANGE   Modify an equioscillation reference using the Remez algorithm.
%   EXCHANGE(XK, H, METHOD, F, P, Q) performs one step of the Remez algorithm
%   for the best rational approximation of the CHEBFUN F of the target function
%   according to the first method (METHOD = 1), i.e. exchanges only one point,
%   or the second method (METHOD = 2), i.e. exchanges all the reference points.
%   XK is a column vector with the reference, H is the levelled error, P is the
%   numerator, and Q is the denominator of the trial rational function P/Q.
%
%   [XK, NORME, E_HANDLE, FLAG] = EXCHANGE(...) returns the modified reference
%   XK, the supremum norm of the error NORME (included as an output argument,
%   since it is readily computed in EXCHANGE and is used later in REMEZ), a
%   function handle E_HANDLE for the error, and a FLAG indicating whether there
%   were at least N+2 alternating extrema of the error to form the next
%   reference (FLAG = 1) or not (FLAG = 0).
%
%   [XK, ...] = EXCHANGE([], 0, METHOD, F, P, Q, N + 2) returns a grid of N + 2
%   points XK where the error F - P/Q alternates in sign (but not necessarily
%   equioscillates). This feature of EXCHANGE is useful to start REMEZ from an
%   initial trial function rather than an initial trial reference.

% Compute extrema of the error.
e_num = diff(f - p) ;
rts = roots(e_num, 'nobreaks');
rr = [f.domain(1) ; rts; f.domain(end)];
if ( length(bands) ==  4 )
    %idx = (rr > dom(2)) & (rr < dom(3))
    %rr(idx) = [];
    rr = sort( [rr; bands(2); bands(3)] );
end

% Function handle output for evaluating the error.
err_handle = @(x) feval(f, x) - feval(p, x);

% Select exchange method.
if ( method == 1 )                           % One-point exchange.
    [ignored, pos] = max(abs(feval(err_handle, rr)));
    pos = pos(1);
else                                           % Full exchange.
    pos = find(abs(err_handle(rr)) >= abs(h)); % Values above levelled error
end

% Add extrema nearest to those which are candidates for exchange to the
% existing exchange set.
[r, m] = sort([rr(pos) ; xk]);
v = ones(Npts, 1);
v(2:2:end) = -1;
er = [feval(err_handle, rr(pos)) ; v*h];
er = er(m);

% Delete repeated points.
repeated = diff(r) == 0;
r(repeated) = [];
er(repeated) = [];

% Determine points and values to be kept for the reference set.
s = r(1);    % Points to be kept.
es = er(1);  % Values to be kept.
for i = 2:length(r)
    if ( (sign(er(i)) == sign(es(end))) && (abs(er(i)) > abs(es(end))) )
        % Given adjacent points with the same sign, keep one with largest value.
        s(end) = r(i);
        es(end) = er(i);
    elseif ( sign(er(i)) ~= sign(es(end)) )
        % Keep points which alternate in sign.
        s = [s ; r(i)];    %#ok<AGROW>
        es = [es ; er(i)]; %#ok<AGROW>
    end
end

% Of the points we kept, choose n + 2 consecutive ones that include the maximum
% of the error.
[norme, index] = max(abs(es));
d = max(index - Npts + 1, 1);
if ( Npts <= length(s) )
    xk = s(d:d+Npts-1);
    flag = 1;
else
    xk = s;
    flag = 0;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions for displaying diagnostic information.

% Function called when opts.plotIter is set.
function doPlotIter(xo, xk, err_handle, dom)

xxk = linspace(dom(1), dom(end), 300);
plot(xo, 0*xo, 'or', 'MarkerSize', 12)   % Old reference.
holdState = ishold;
hold on
plot(xk, 0*xk, '*k', 'MarkerSize', 12)   % New reference.
plot(xxk, err_handle(xxk))               % Error function.
if ( ~holdState )                        % Return to previous hold state.
    hold off
end
xlim(dom)
legend('Current Ref.', 'Next Ref.', 'Error')
drawnow

end

% Function called when opts.displayIter is set.
function doDisplayIter(iter, err, h, delta, normf, diffx)

disp([num2str(iter), '        ', num2str(err, '%5.4e'), '        ', ...
    num2str(abs(h), '%5.4e'), '        ', ...
    num2str(delta/normf, '%5.4e'), '        ', num2str(diffx, '%5.4e')])

end