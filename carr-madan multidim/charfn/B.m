function B = B(tau, a, b, kappa, sigma, theta)
    [d, gamma, ~, ~] = get_d_and_gamma(kappa, sigma, a, b, tau);
    B = kappa * theta / sigma ^ 2 * (kappa - d) * tau + ....
        2 * kappa * theta / sigma ^ 2 * log(2 * d / gamma);
end