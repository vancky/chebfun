function varargout = chebfir(n, freqs, f, varargin)
%CHEBFIR   Chebfun filter desigin for FIR filters.
%
% See also CF.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.
if ( ~isa(f, 'chebfun') )
    f = chebfun(@(x), e
end

varargout = remez(f, n, freqs, varargin);
end