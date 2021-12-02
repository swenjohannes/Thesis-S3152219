function B2 = B2(w1, w2, kappa, rho, eta, tau)

eta2 = eta ^ 2;
beta = kappa - 1i * rho * eta * w1;
D = sqrt(beta ^ 2 + eta2 * w1 * (1i + w1));

h = (beta - D - 1i * w2 * eta2) / (beta + D - 1i * w2 * eta2);
hemdt = h  * exp (-D * tau); %To ease up calculations

%Return values
B2 = (beta - D - (beta + D) * hemdt) / (eta2 * (1 - hemdt));
end

