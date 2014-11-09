% compute roots of a trigtech
function r = jroots2(c)

% find intervals with roots
root_ints = zero_ints(c);

% compute roots using newton's method
rts = newton_roots(c,root_ints);

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
minc = 10;

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

% newton steps
for ii=1:5
   stp = c.feval(rts)./dc.feval(rts);
   norm(stp,inf)
   rts = rts - stp;    
end

end
