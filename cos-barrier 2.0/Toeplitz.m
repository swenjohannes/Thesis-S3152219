%Toeplitz multiplication via FFT
c = [3 4 5 6];
r = [3 0 9 2];

M = toeplitz(c, r);
u = [2 4 6 8];
(M * u')'

v = [c r(1) flip(r(2:4))];
us = [u repelem(0, 4)];

res = ifft(fft(v) .* fft(us));
res(1:4)

%Now complex:
ci = c + i * [5 0 1 2];
ri = r + i * [5 1 7 8];

Mi = toeplitz(ci, ri);
ui = u + i * [1 2 3 4];
(Mi * ui')'

vi = [ci ri(1) flip(ri(2:4))];
usi = [ui repelem(0, 4)];
resi = ifft(fft(vi) .* fft(usi));
resi(1:4)
