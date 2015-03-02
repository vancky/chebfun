function A = uminus(A)
    %UMINUS   CHEBOP unary minus.
    A.op = @(varargin) - A.op(varargin{:});
end
