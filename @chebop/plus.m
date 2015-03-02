function out = plus(A,B)
    %PLUS   CHEBOP plus.
    if ( isa(A, 'chebop') && isa(B, 'chebop') )
        if ( length(A.inDims) == length(B.inDims) ...
             && all(A.inDims == B.inDims) )
            % Make sure sizes match.
            out = chebop(@(varargin) A.op(varargin{:}) + B.op(varargin{:}));
            out.inDims = A.inDims;
        else
            error('CHEBFUN:CHEBOP:plus:sizeMismatch', ...
                'Operator dimensions must agree.')    
        end

    elseif ( isa(A, 'chebfun') || isa(A, 'double') )
        out = chebop(@(varargin) B.op(varargin{:})+A);

    elseif ( isa(B, 'chebfun') || isa(B, 'double') )
        out = chebop(@(varargin) A.op(varargin{:})+B);

    else
        if ( ~isa(A, 'chebop') )
            objType = class(A);
        else
            objType = class(B);
        end
        error('CHEBFUN:CHEBOP:plus:notSupported', ...
            ['CHEBOP/PLUS not supported for object of type ' ...
            objType '.'])
    end
end
