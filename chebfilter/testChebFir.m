n = 3;
Fp = .4;
Fs = .6;
Fc = (Fp + Fs)/2;
Ft = (Fp - Fs);
tic
[h, dm] = firpm(2*n, [0 Fp Fs 1], [1 1 0 0 ] );
toc
Hideal = @(x) 1.*(x>=0) - 1.*(x>.5);
tic
[p, dc]= chebfir( n, [0, .4, .6, 1], Hideal );
toc
h = h(:);
dom = [-1, 1];
a = dom(1);
b = dom(2);
T = b - a;
H = chebfun(h, dom, 'periodic', 'coeffs');
dc - dm
norm(p.coeffs - h, inf)
d = dm;
%%
% Plotting
subplot(2, 1, 1)
plot(H)
hold on
plot(p, '.r')
w = chebfun(@(w) w, dom);
u = heaviside(w+T/4)-heaviside(w-T/4);
plot(u)
plot(T/2*[-Fp, Fp], 1 + d*[1, 1], '--g')
plot(T/2*[-Fp, Fp], 1 - d*[1, 1], '--g')
plot(T/2*[-Fp, -Fp, NaN, Fp, Fp], [1-d, 1+d, NaN, 1-d, 1+d], '--g' )

plot(T/2*[-1, -Fs, NaN, Fs, 1], d*[1, 1, NaN, 1, 1], '--r')
plot(T/2*[-1, -Fs, NaN, Fs, 1], -d*[1, 1, NaN, 1, 1], '--r')
plot(T/2*[-Fs, -Fs, NaN, Fs, Fs], [-d, d, NaN, -d, d], '--r' )
hold off
%Compare with CF
subplot(2, 1, 2)
% Random points in the pass band:
ff = sort(T/2*(-Fp + 2*Fp*rand(1000,1)));
plot(ff, p(ff)-H(ff))
%hold on
% append random points in the stop band:
ff = sort(T/2*(Fs + (1-Fs)*rand(500,1))); 
plot(ff, p(ff)-H(ff));
ff = sort(T/2*(-1 + (1-Fs)*rand(500,1)));
plot(ff, p(ff)-H(ff));
plot(p - H)
hold off
[p.coeffs , H.coeffs]
%%
% MATLAB style
% [H, W] = freqz(h, 1, 1000);
% subplot(2, 1, 1)
% plot(W, abs(H));
% subplot(2, 1, 2);
% plot(W, angle(H));