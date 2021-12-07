function F = get_F(x, v, kappa, rho, eta, theta, r, tau, a1, b1, a2, b2, N)

pbma1 = pi / (b1 - a1);
pbma2 = pi / (b2 - a2);

[Fp, Fm] = deal(zeros(N));
for k1 = 0:(N - 1)
    for k2 = 0:(N - 1)
        Fp(k1 + 1, k2 + 1) = phi(x, v, k1 * pbma1, k2 * pbma2, kappa, rho, eta, theta, r, tau) ...
                                .* exp(-1i * k1 * pbma1 * a1  -  1i * k2 * pbma2 * a2);
        Fm(k1 + 1, k2 + 1) = phi(x, v, k1 * pbma1, -k2 * pbma2, kappa, rho, eta, theta, r, tau) ...
                                .* exp(-1i * k1 * pbma1 * a1 +  1i * k2 * pbma2* a2);
    end
end
F = real(Fp + Fm);
end
