Fp = .445;
Fs = .542;
Fc = (Fs+Fp)/2;
n = 24;
freqs = [0, Fp, Fs, 1];
f = chebfun(@(x) heaviside(Fc-x), 'splitting', 'on');
Hc = firpm(n, freqs, f);
[h, dm, res] = firpm(n, freqs, [1 1 0 0 ] );
h = h(:);
%%
% We now make a periodic chebfun from the filter
% coefficients computed:
if ( rem(n, 2) == 0 )
    if( norm(h([1, end]), inf) < eps )
        h = h(2:end-1);
    end
end
Hm = chebfun(h, 'coeffs', 'trig');
%%
subplot(3, 1, 1)
plot(Hc, 'r')
hold on
plot(Hm)
hold off
subplot(3, 1, 2)
ec = filterError(f, Hc, Fp, Fs);
em = filterError(f, Hm, Fp, Fs);
plot(ec, 'r')
hold on
plot(em, 'k')
hold off
subplot(3, 1, 3)
fprintf('          %s                %s\n', 'Chebfun' , 'Matlab');
s = [minandmax(ec),  minandmax(em)];
disp(s)