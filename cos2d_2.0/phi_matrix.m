
phi_ = @(w1, w2) phi(x, v, w1 * pbma1, w2 * pbma2, kappa, rho, eta, theta, r, tau);

phi_p = zeros(N);
phi_m = zeros(N);
for k1 = 0:(N - 1)
    for k2 = 0:(N - 1)
        phi_p(k1 + 1, k2 + 1) = phi_(k1, k2);
        phi_m(k1 + 1, k2 + 1) = phi_(k1, -k2);
    end
end