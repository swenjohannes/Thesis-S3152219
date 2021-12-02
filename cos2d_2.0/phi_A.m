function phi_A = phi_A(w1, w2, kappa, rho, eta, theta, r, tau)

eta2 = eta ^ 2;
beta = kappa - 1i * rho * eta * w1;
D = sqrt(beta ^ 2 + eta2 * w1 * (1i + w1));

h = (beta - D - 1i * w2 * eta2) / (beta + D - 1i * w2 * eta2);
hemdt = h  * exp (-D * tau); %To ease up calculations

A  = 1i * w1 * r * tau + ... %Compute A
    kappa * theta / eta2 * ((beta - D) * tau - 2 * log((hemdt - 1) / (h - 1)));
phi_A = exp(A);
end