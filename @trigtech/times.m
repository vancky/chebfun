function f = times(f, g, varargin)
%.*   TRIGTECH multiplication.
%   F.*G multiplies TRIGTECH objects F and G or a TRIGTECH by a scalar if either
%   F or G is a scalar.
%
%   If F is an array-valued TRIGTECH, then F.*C is supported if C is a row
%   vector of doubles with the same number of columns as F.
%
% See also MTIMES, RDIVIDE.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% TRIGTECH * [] = []:
if ( isempty(f) || isempty(g) )
    f = []; 
    return
end

if ( ~isa(f, 'trigtech') )      % Ensure F is a TRIGTECH.
    
    f = times(g, f, varargin{:});
    return
    
elseif ( isa(g, 'double') )     % TRIGTECH .* double.
    
    % Do the multiplication:
    if ( size(g, 2) > 1 )
        f.coeffs = bsxfun(@times, f.coeffs, g);
        f.values = f.coeffs2vals(f.coeffs);
        f.vscale = f.vscale.*abs(g);
    else
        f.values = f.values*g;
        f.coeffs = f.coeffs*g;
        f.vscale = f.vscale*abs(g);
    end
    f.epslevel = f.epslevel + eps(g);
    f.isReal = f.isReal & isreal(g);
    return

elseif ( ~isa(f, 'trigtech') || ~isa(g, 'trigtech') )
    % Don't know how to do the operation.
    
    error('CHEBFUN:TRIGTECH:times:typeMismatch', ...
        ['Incompatible operation between objects.\n', ...
         'Make sure functions are of the same type.']);
    
elseif ( size(f.values, 1) == 1 )
    % If we have (constant TRIGTECH).*TRIGTECH, reverse the order and call TIMES
    % again:
    f = times(g, f.coeffs);
    f.epslevel = max(f.epslevel, g.epslevel);
    return
    
elseif ( size(g.values, 1) == 1)
    % If we have TRIGTECH.*(constant TRIGTECH), convert the (constant TRIGTECH)
    % to a scalar and call TIMES again:
    f = times(f, g.coeffs); 
    f.epslevel = max(f.epslevel, g.epslevel);
    return
end

% Do muliplication in coefficient space:
[f.coeffs, pos] = coeff_times_main(f.coeffs, g.coeffs);
f.values = f.coeffs2vals(f.coeffs);

% Update vscale, epslevel, and ishappy:
vscale = max(abs(f.values), [], 1);

% Avoid NaNs:
tmpVscale = vscale;
tmpVscale(vscale == 0) = 1;
f.vscale(f.vscale == 0) = 1;
g.vscale(g.vscale == 0) = 1;
 
% See CHEBTECH CLASSDEF file for documentation on this:
epslevelBound = (f.epslevel + g.epslevel) .* (f.vscale.*g.vscale./tmpVscale);
f.epslevel = updateEpslevel(f, epslevelBound);
f.vscale  = vscale;
f.ishappy = f.ishappy && g.ishappy;

% Simplify.
f = simplify(f);

if ( pos )
    % Here we know that the product of F and G should be positive. However,
    % SIMPLIFY may have destroyed this property, so we enforce it.
    f.values = abs(f.values); 
    f.coeffs = f.vals2coeffs(f.values);
    f.isReal = true(1, size(f.coeffs, 2));
else
    f.isReal = f.isReal & g.isReal;
end

% If f and g are real then make the result real.
f.values(:,f.isReal) = real(f.values(:,f.isReal));

end

function [coeffs, pos] = coeff_times_main(f, g)

% Get the size of each TRIGTECH:
[fn, fm] = size(f);
[gn, gm] = size(g);

% Prolong:
f((fn+1):(fn+gn-1),:) = 0;
g((gn+1):(fn+gn-1),:) = 0;

% Check dimensions:
if ( fm ~= gm )
    if ( fm == 1 )
        % Allow [Inf x 1] .* [Inf x m].
        f = repmat(f, 1, gm);
    elseif ( gm == 1 )
        % Allow [Inf x m] .* [Inf x 1].
        g = repmat(g, 1, fm);
    else
        error('CHEBFUN:TRIGTECH:times:dim2', ...
            'Inner matrix dimensions must agree.');
    end
end

% Check for two cases where the output is known in advance to be positive,
% namely F == conj(G) or F == G and isreal(F).
pos = false;

% Multiply coefficients:
if ( isequal(f, g) )
   coeffs = coeff_times(f, g);          
   if ( isreal(f) )
       pos = true; 
   end
elseif ( isequal(conj(f), g) )
   coeffs = coeff_times(f, g);
   pos = true;
else
   if fn > gn
       coeffs = coeff_times(f, g);
   else
       coeffs = coeff_times(g, f);
   end
end

end

function hc = coeff_times(fc, gc)
%COEFF_TIMES(FC, GC)   Multiplication in coefficient space.
%   HC = COEFF_TIMES(FC, GC) returns the vector of Fourier coefficients, HC,
%   resulting from the multiplication of two functions with FC and GC
%   coefficients. The vectors have already been prolonged.

N = length(fc);
M = length(gc);
gc = gc(1:M,:);
multmat = trigspec.multmat(N+M-1,fc);
multmat = multmat(:,1:M);
hc = multmat*gc;    

%hc = ifft(fft(fc).*fft(gc));

% mn = length(fc);
% t = [2*fc(1,:) ; fc(2:end,:)];                    % Toeplitz vector.
% x = [2*gc(1,:) ; gc(2:end,:)];                    % Embed in Circulant.
% xprime = fft([x ; x(end:-1:2,:)]);                % FFT for Circulant mult.
% aprime = fft([t ; t(end:-1:2,:)]);
% Tfg = ifft(aprime.*xprime);                   % Diag in function space.
% hc = .25*[Tfg(1,:); Tfg(2:end,:) + Tfg(end:-1:2,:)];% Extract out result.
% hc = hc(1:mn,:);  

end
