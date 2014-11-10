% compute roots of a trigtech
function r = jroots2(c)

% find intervals with roots
root_ints = zero_ints(c);

% compute roots using newton's method
rts = newton_roots(c,root_ints);
% rts = secant_roots(c,root_ints);

% X = linspace(-1,1,401)';
% 
% hold on;
% % x-axis
% plot(X,0*X,'k');
% plot(root_ints,0*root_ints,'ro');
% hold off;

r = rts;

end

function root_ints = zero_ints(c)

% max density increase
minc = 100;

% set maximum power of 2
mpow = 20;

% get coeffs
coeffs = c.coeffs;

% get degree
N = length(coeffs);

% set new length as a power of 2
m = min(mpow,ceil(log2(minc*N)));
% m = mpow;
if mod(N,2) == 0
    M = 2^m;
else
    M = 2^m+1;
end 

% set new coeffs
if mod(N,2) == 0
    new_cfs = [zeros((M-N)/2,1);coeffs(1)/2;coeffs(2:end);coeffs(1)'/2;zeros((M-N)/2-1,1)];
else
    new_cfs = [zeros((M-N)/2,1);coeffs;zeros((M-N)/2,1)];
end

% compute endpoints
intervals = linspace(-1,1,M+1);
intervals = [intervals(1:end-1)',intervals(2:end)'];

% compute function values
vals = trigtech.coeffs2vals(new_cfs);
vals = [vals;vals(1)];

% compute signs
val_signs = sign(vals);

% find intervals that contain roots
root_ints = intervals((abs(diff(val_signs))==2),:);

% set new_numints
new_num_ints = size(root_ints,1);

end


function rts = newton_roots(c,root_ints)

% initial guesses
rts = sum(root_ints,2)/2;

% df
dc = diff(c);

STPS = [];
FS = [];
DFS = [];
% newton steps
for ii=1:5
   stp = c.feval(rts)./dc.feval(rts);
%    STPS = [STPS,stp];
%    FS = [FS,c.feval(rts)];
%    DFS = [DFS,dc.feval(rts)];
   nrm = norm(stp,inf)
   if nrm < eps
       break;
   end
   rts = rts - stp;    
end
% STPS
% FS
% DFS

end

function rts = secant_roots(c,root_ints)

% initial guesses
old_rts = sum(root_ints,2)/2;
rts = sum([root_ints(:,1),old_rts],2)/2;

% df
dc = diff(c);

% secant steps
for ii=1:20
   stp = c.feval(rts).*(rts-old_rts)./(c.feval(rts)-c.feval(old_rts));
   stp(isnan(stp)) = 0;
   stp(isinf(stp)) = 0;
   nrm = norm(stp,inf)
   dnrm = norm(dc.feval(rts),inf);
   if nrm < eps*dnrm
       break;
   end
   old_rts = rts;
   rts = rts - stp;
end

end
