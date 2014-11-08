% compute roots of a trigtech
function r = jroots2(c)

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
if mod(N,2) == 0
    M = 2^m;
else
    M = 2^m+1;
end 

% set new coeffs
if mod(N,2) == 0
    newcfs = [zeros((M-N)/2,1);coeffs;coeffs(1)';zeros((M-N)/2-1,1)];
else
    newcfs = [zeros((M-N)/2,1);coeffs;zeros((M-N)/2,1)];
end

% compute function values
vals = trigtech.coeffs2vals(newcfs);
vals = [vals;vals(1)];

% compute signs
valsigns = sign(vals);
intervals = [valsigns,valsigns([2:M+1,1])];

% find intervals that contain roots
rootinds = (abs(sum(intervals,2))==0);

r = valsigns;


end
