function sol = mldivide(A, rhs)
%MLDIVIDE   Solve a continuous system of nonlinear equations.

sizesA = blockSizes(A);
sizesRHS = blockSizes(rhs);

% NOTE: This code assumes that the chebmatrices A and RHS are well-formed.
% The chebmatrix class itself is responsible for making sure concatenations
% and so on don't produce malformed chebmatrices.

% TODO: Check to make sure the RHS and the operator A are of consistent
% dimensions.

% TODO: Check that the domains of A and RHS are consistent.
dom = rhs.domain;

% % Convert chebfun rows to functionalBlocks.
% s = @(f) functionalBlock.inner(f);
% for i = 1:size(A, 1),
%     for j = 1:size(A, 2)
%         u = A.blocks{i,j};
%         if ( isa(u, 'chebfun') )
%             if ( u.isTransposed )
%                 u = functionalBlock.inner(u.', dom);
%             end
%         end
%     end
% end

% First treat the operator part.
isOp = cellfun(@(f) isa(f, 'chebop'), A.blocks);
if ( sum(isOp(:)) == 1 )
    % There is only one operator.
    N = A.blocks{isOp};
    op = N.op;

    % TODO: Take care of column chebfuns in the same row as N.

    % Take care of boundary conditions. We know they are just a set of rows,
    % so loop over them. We need to modify the RHS accordingly.
    isBC = ~isOp;
    bc = '@(t,u)[';
    for k = find(isBC).'
        bc = [bc ...
          'sum((A.blocks{' num2str(k) '}.'').*u)-rhs.blocks{' num2str(k) '};'];
    end
    bc = str2func([bc ']']);
    % rhs.blocks{isBC} = [];  % Why doesn't this work?
    rhs.blocks = subsasgn(rhs.blocks, ...
        struct('type', '()', 'subs', {{find(isBC)}}), []);
else
    % Not yet.
    error('CHEBFUN:CHEBMATRIX:mldivide:unsupported', ...
        'Solving systems with multiple chebops is not yet supported.')
end

B = chebop(op, dom, bc);
sol = B \ rhs;

end
