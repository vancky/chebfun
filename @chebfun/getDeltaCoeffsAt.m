function DV = getDeltaCoeffsAt(F, loc, N)
%GETDELTACOEFFSAT   Delta and delta derivative coeffs at a location.
%   GETDELTACOEFFSAT(F,LOC,N) returns a vector of the coefficients of the
%   first N-1 derivatives of delta functions at the point LOC in F.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% TODO: Check that LOC is in the domain of F.

if ( N < 0 )
    DV = [];
    return
end

DV = zeros(N, 1);

% I'm pretty sure this is all going to have to be rewritten to not
% break encapsulation, but.
for k = 1:length(F.funs)
    fun = F.funs{k};

    if ( strcmpi(class(fun), 'deltafun') )
        % We have to check that this fun is a deltafun, since otherwise
        % it doesn't have the deltaLoc property.
        locs = F.funs{1}.deltaLoc;

        if ( any(find(abs(locs - loc) < eps)) )
            % Find the relevant delta magnitudes.
            col = find(abs(locs - loc) < eps);
            mag = fun.deltaMag(:,col);

            if ( size(mag,1) > length(DV) )
                % We don't know which is larger.
                DV = mag(1:N);
                return
            else
                DV(1:length(mag)) = mag;
            end

        end

    else
        % The current fun is not a deltafun.
        continue
    end
end

end
