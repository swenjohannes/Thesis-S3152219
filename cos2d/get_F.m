function F = get_F(x, v, kappa, rho, eta, theta, r, tau, a, b, N)
bma = b - a;
pbma = pi / bma;

phi_ = @(w1, w2) phi(x, v, w1 * pbma, w2 * pbma, kappa, rho, eta, theta, r, tau);

phi_p = zeros(N);
phi_m = zeros(N);
for k1 = 0:(N - 1)
    for k2 = 0:(N - 1)
        phi_p(k1 + 1, k2 + 1) = phi_(k1, k2);
        phi_m(k1 + 1, k2 + 1) = phi_(k1, -k2);
    end
end

%To efficiently compute the matrices  use k1' + k2
k1 = 0:(N - 1);
k2 = 0:(N - 1);
Fp = phi_p .* exp(-1i * (k1' + k2) * a * pbma);
Fm = phi_m .*exp(-1i * (k1' - k2) * a * pbma);
F = real(Fp + Fm);
end