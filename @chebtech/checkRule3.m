function [ishappy, epslevel, cutoff] = checkRule3(f, values, pref)
%CLASSICCHECK   Attempt to trim trailing Chebyshev coefficients in a CHEBTECH.
%   [ISHAPPY, EPSLEVEL, CUTOFF] = CLASSICCHECK(F, VALUES) returns an estimated
%   location, the CUTOFF, at which the CHEBTECH F could be truncated to
%   maintain an accuracy of EPSLEVEL relative to F.VSCALE and F.HSCALE. ISHAPPY
%   is TRUE if the representation is "happy" in the sense described further
%   below and FALSE otherwise.  If ISHAPPY is FALSE, EPSLEVEL returns an
%   estimate of the accuracy achieved.


% CALCULATE THE NORMALIZED ENVELOPE SEQUENCE
c = f.coeffs;
b = abs(c);
m = b(end)*ones(size(b));
for j = length(b)-1:-1:1
  m(j) = max(b(j),m(j+1));
end
a = m/m(1);

% call rule 3
[epslevel, cutoff] = rule3(a);

% set ishappy
ishappy = false;
if cutoff < length(c)
    ishappy = true;
end

end


function [epslevel, cutoff] = rule3(a)

%%
% A chopping rule of "standard" type, that is, with an input
% tolerance (currently hardwired) that is applied with some
% flexibility.  Nick Trefethen, 24 November 2014.

%%
% INPUT: a sequence a with length >= 17, starting with a(1) = 1,
% monotonically nonincreasing.
%
% OUTPUT: a positive integer CUTOFF.  If CUTOFF = length(a), then rule 3
% is "not happy": a satisfactory chopping point has not been found.
% If CUTOFF < length(a), the rule is "happy" and CUTOFF represents the
% last index in the sequence that should be retained.

%%
% This rule works from a tolerance TOL, which by default is
% machine epsilon (2.2e-16).  The series will never be chopped
% unless it falls at least below TOL^(2/3).  It will always be
% chopped if it has a long enough segment below TOL.

  n = length(a);
  cutoff = n;
  epslevel = a(cutoff);
  if n < 17, return, end
  tol = 2^(-52);

%%
% Step 1: Scan the sequence for a "plateau point K", the first
% point j, if any, that is followed by a plateau.  A plateau is a
% stretch of coefficients a(j),...,a(j2), j2 = round(1.25*j+5) <= n,
% with the property that a(j2)/a(j) > r.  The number r ranges
% from r = 0 if a(j) = TOL to r = 1 if a(j) = TOL^(2/3).
% (Thus a plateau with a(j) ~ TOL^(2/3) has to be perfectly flat
% to count, whereas for a(j) ~ TOL it doesn't have to be flat at all.)
% If a plateau point K is found, then we know we are going to chop the
% series, but the precise chopping point J >= K-1 is not yet determined.

  for j = 1:n
    j2 = round(1.25*j+5); 
    if j2 > n, return, end         % there is no plateau: exit
    a1 = a(j);
    a2 = a(j2);
    r = 3*(1-log(a1)/log(tol));
    plateau = (a1 == 0) | (a2/a1 > r);
    if plateau, K = j; break, end
  end

%%
% Step 2: Chop the sequence at the point of the plateau where
% the coefficient sequence a, plus a linear function included
% to bias the result towards the left end, is minimal.
  
  if a(K) == 0
      cutoff = K-1;
  else
      aa = log10(a(K:j2)); aa = aa(:);
      aa = aa + (1+.25*log10(a(K-1)/a(K)))*linspace(0,1,j2+1-K)';
      [ignore,d] = min(aa);
      cutoff = K-2+d;
  end

  epslevel = a(cutoff+1);
end

