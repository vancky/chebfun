function xk = getInitialReference(freqs, n, N)
% n is the order of the filter and N is the number of points in the exchange
% set. N = n.

% If the above procedure failed to produce a reference with enough 
% equioscillation points, just use the Chebyshev points.
L = length(freqs);
xk = [];

% Special case:
if ( L == 4 )
    xk = [linspace(freqs(1), freqs(2), ceil((N+2)/2)).'; ...
        linspace(freqs(3), freqs(4), floor((N+2)/2)).'; ];
else
    % L must be even
    for k = 1:L/2
        a = freqs(2*k-1);
        b = freqs(2*k);
        nk = round((N+2)*(b-a)/pi);
        xk = [xk; linspace(a, b, nk).'];
    end
end

if ( length(xk) > N + 2 )
    % TODO
    error( 'N+2 length not achieved.');
end

if ( length(xk) < N + 2 )
    % TODO
    error( 'N+2 length not acheived.');
end
    
end
