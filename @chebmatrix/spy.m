function spy(A)
%SPY    Visualize a chebmatrix.
%   If A is a chebmatrix, SPY(A) creates a picture of the nonzero pattern of the
%   default discretization of A. Block boundaries are indicated by gray
%   lines. 
%   
%   SPY(A,DIM) uses the dimension vector DIM to create the picture.
%
%   See also CHEBMATRIX, CHEBMATRIX.MATRIX, CHEBOPPREF.

%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.

% Select a reasonable default discretization.
dom = A.domain;
if length(dom)==2   % no breakpoints
    dim = 10;
else
    dim = repmat(6,[1 length(dom)-1]);
end

% Discretize and do a regular spy plot.
data = matrix(A,dim,dom);
spy(data)
s =  sprintf('%i,',dim);    % list of sizes
s = [ 'discretization = [', s(1:end-1), ']' ];
xlabel(s)

% Override hold state.
hs = ishold;
hold on

% Find all block sizes, substituting in the discretization size for Inf.
[m,n] = blockSizes(A);
m(isinf(m)) = sum(dim);  
n(isinf(n)) = sum(dim);

% Draw horizontal block boundaries.
csrow = cumsum(m(:,1)');
rowdiv = csrow(1:end-1) + 1/2;   % insert boundary after each block
colmax = sum(n(1,:));
plot([0;colmax+1],[rowdiv;rowdiv],'color',[.6 .6 .6])

% Draw vertical block boundaries.
cscol = cumsum(n(1,:));
coldiv = cscol(1:end-1) + 1/2; 
rowmax = sum(m(:,1));
plot([coldiv;coldiv],[0;rowmax+1],'color',[.6 .6 .6])

% Clean up hold state.
if ( ~hs )
    hold off
end

end