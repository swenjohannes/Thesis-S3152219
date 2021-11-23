function A = A(tau, a, b, kappa, sigma)

% obtain the values for d, gamma and c1/c2   
[d, gamma, c1, c2] = get_d_and_gamma(kappa, sigma, a, b, tau);

numerator = d * a * c1 - c2 * (2 * b + kappa * a);
A = numerator / gamma; %Compute A

end