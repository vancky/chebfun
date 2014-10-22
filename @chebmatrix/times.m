function C = times(A, B)
%.*   Elementwise products with CHEBMATRICES.
%
% See also TIMES.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Scalar operand is scalar*identity. But it's faster to interpret it as a
% special case.

if ( isnumeric(A) && (numel(A)==1) )
    C = scalarTimes(B, A);

elseif ( isnumeric(B) && (numel(B)==1) )
    C = scalarTimes(A, B);

else
    C = elementwiseTimes(A, B);
end

end


function C = scalarTimes(A, z)
%SCALARTIMES   Multiplying blocks in a CHEBMATRIX with a scalar.
[m, n] = size(A);
C = cell(m, n);
for i = 1:m
    for j = 1:n
        C{i,j} = z * A.blocks{i,j};
    end
end
C = chebmatrix(C);
end


function C = elementwiseTimes(A, B)
%ELEMENTWISETIMES   Multiply CHEBMATRIX .* other (or the other way around).

% All the row-entries of A must be of the same type, so we HORZCAT. This is much
% much much faster than converting a matrix to a CHEBMATRIX and looping over all
% the block entries.
m = size(A,1);
C = cell(m, size(B,2));
for k = 1:m
    tmp = horzcat(A.blocks{k,:}); % Horzcat
    if ( isa(tmp, 'chebfun') )    % Use array-valued if possible
        tmp = quasi2cheb(tmp);    % Since we're adding, all breaks
    end                           % will eventually be the same.
    if ( isvector(B) )
        C(k) = tmp*B;             % Avoid num2cell call.
    else
        C(k,:) = num2cell(tmp*B);
    end
end
% Convert cell to CHEBMATRIX for ouput:
C = chebmatrix(C);
end
