function varargout = chebfir(n, freqs, f, varargin)
%CHEBFIR   
%
%   [1] Pachon, R. and Trefethen, L. N.  "Barycentric-Remez algorithms for best
%   polynomial approximation in the chebfun system", BIT Numerical Mathematics,
%   49:721-742, 2009.
%
%   [2] Pachon, R.  "Algorithms for Polynomial and Rational Approximation".
%   D. Phil. Thesis, University of Oxford, 2010 (Chapter 6).
%
% See also CF.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% To be used by the algorithm as the denominator:
q = chebfun(1);

[N, freqs, f, opts] = parseInputs(n, freqs, f, varargin{:});
normf = norm(f);
% Initial values for some parameters.
iter = 0;       % Iteration count.
delta = normf;  % Value for stopping criterion.
deltamin = inf; % Minimum error encountered.
diffx = 1;      % Maximum correction to trial reference.

% Compute an initial reference set to start the algorithm.
vk = getInitialReference(freqs, n, N);

% Transform the problem to that of polynomial case:
xk = cos(pi*vk);
fx = chebfun(@(x) f(1/pi*acos(x)), 'splitting', 'on');
xo = xk;

% Print header for text output display if requested.
if ( opts.displayIter )
    disp('It.     Max(|Error|)       |ErrorRef|      Delta ErrorRef      Delta Ref')
end

% Run the main algorithm.
while ( (delta/normf > opts.tol) && (iter < opts.maxIter) && (diffx > 0) )
    fk = feval(fx, xk);     % Evaluate on the exchange set.
   
    % Barycentric weights for exchange set.
    w = baryWeights(xk);
        
    [p, h] = computeTrialFunctionPolynomial(fk, xk, w, n, N);
    
    % Perturb exactly-zero values of the levelled error.
    if ( h == 0 )
        h = 1e-19;
    end

    % Update the exchange set using the Remez algorithm with full exchange.
    [xk, err, err_handle] = exchange(xk, h, 2, fx, p, q, N + 2, freqs);

    % If overshoot, recompute with one-point exchange.
    if ( err/normf > 1e5 )
        [xk, err, err_handle] = exchange(xo, h, 1, fx, p, q, N + 2, freqs);
    end

    % Update max. correction to trial reference and stopping criterion.
    diffx = max(abs(xo - xk));
    delta = err - abs(h);

    % Store approximation with minimum norm.
    if ( delta < deltamin )
        pmin = p;
        if ( n > 0 )
            qmin = q;
        end

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
    warning('CHEBFUN:CHEBFUN:remez:convergence', ...
        ['Remez algorithm did not converge after ', num2str(iter), ...
         ' iterations to the tolerance ', num2str(opts.tol), '.']);
end

% Form the outputs.
status.delta = delta/normf;
status.iter = iter;
status.diffx = diffx;
status.xk = xk;

%p = simplify(p);
c = chebcoeffs(p);
% Convert the cosine coefficients to complex exponentials:
a = [1/2*flipud(c(2:end)); c(1); 1/2*c(2:end)];
p = chebfun( a, [-1, 1], 'periodic', 'coeffs');
varargout = {p, err, status};

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input parsing.

function [N, freqs, f, opts] = parseInputs(n, freqs, f, varargin)

if ( n < 1 || n ~= round(n) || ~isreal(n) )
    error('CHEBFUN:CHEBFUN:chebfir:badInput', ...
            'n must be a positive integer.')
end

if ( min(size(freqs)) > 1 || any(diff(freqs) <= 0) || any((freqs < 0) | (freqs >  pi)) )
    error('CHEBFUN:CHEBFUN:chebfir:badInput', ...
            'Frequencies must be an ascending vector with elements in [0, pi]')
end

if ( isa(f, 'chebfun') )
    if ( any(f.domain ~= [0, 1]) )
        error('CHEBFUN:CHEBFUN:chebfir:badInput', ...
            'Domain of desired response should be [-pi, pi].')
    end
elseif ( isa(f, 'function_handle') )
    f = chebfun(f, [0, 1], 'splitting', 'on');
else
    error('CHEBFUN:CHEBFUN:chebfir:badInput', ...
            'Desired response should be a chebfun or a function handle.');
end
    
% The exchange uses N + 2 points, so for the periodic case, 
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
        error('CHEBFUN:CHEBFUN:chebfir:badInput', ...
            'Unrecognized sequence of input parameters.')
    end
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
