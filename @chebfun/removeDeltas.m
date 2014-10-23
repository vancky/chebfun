function F = removeDeltas(F)
%REMOVEDELTAS   Remove the deltas from a chebfun.
%   REMOVEDELTAS(F) returns a chebfun identical to F with all deltas removed.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

for k = 1:length(F.funs)
    if ( isa(F.funs{k}, 'deltafun') )
        F.funs{k} = F.funs{k}.funPart;
    end
end

F.pointValues = chebfun.getValuesAtBreakpoints(F.funs);

end
