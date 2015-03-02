function out = minus(A,B)
    %MINUS   CHEBOP minus.
    if ( (isa(A, 'chebop') && isa(B, 'chebop')) ...
            || isa(A, 'chebfun') || isa(A, 'double') ...
            || isa(B, 'chebfun') || isa(B, 'double') )

        out = A + (-B);

    else
        if ( ~isa(A, 'chebop') )
            objType = class(A);
        else
            objType = class(B);
        end
        error('CHEBFUN:CHEBOP:minus:notSupported', ...
            ['CHEBOP/MINUS not supported for object of type ' ...
            objType '.'])
    end
end
