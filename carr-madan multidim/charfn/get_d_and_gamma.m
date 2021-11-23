function [d, gamma, c1, c2] = get_d_and_gamma(kappa, sigma, a, b, tau)
    d = sqrt(kappa ^ 2 + 2 * sigma ^ 2 * b);
    c1 = 1 + exp(-d * tau);
    c2 = 1 - exp(-d * tau);
    gamma = d * c1  + (kappa - sigma ^ 2 * a) * c2;
end
