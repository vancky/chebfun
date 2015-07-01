function [u,res] = minres(L,f,tol,maxits)
% Function MINRES

% Check tol
if isempty(tol)
  p = cheboppref;
  tol = p.errTol;
end

% check maxits
if isempty(maxits)
  maxits = 50;
end

% Initialize solution
u = 0*f;

% Initialize Lanczos vectors
[Q,res] = qr(f);
if res == 0, return, end

% Initialize H 
H = [];
D = [];
S = [];

% Iteration loop
for ii = 1:maxits
  
  % Apply operator
  v = L(Q(:,end));

  % orthogonalize
  d1 = Q'*v; v = v - Q*d1;

  % reorthogonalize
  d2 = Q'*v; v = v - Q*d2;

  % normalize
  s1 = norm(v); 
  if s1 == 0
    warning('CHEBFUN/MINRES: operation terminated before convergence')
    break
  end

  % update H
  H(1:ii+1,ii) = [d1+d2;s1];

  % update S and D
  D = [D;H(ii,ii)];
  S = [S;H(ii+1,ii)];

  % construct T
  T = spdiags([S,D,[0;S(1:end-1)]],-1:1,ii+1,ii);

  % solve least squares problem
  u = Q*(T\[res(1);zeros(ii,1)]);

  % compute residual
  res = [res;norm(f-L(u))];
  if res(end) < tol*res(1)
    break
  end

  % update Q
  Q = [Q,v/s1];

end
