function B2 = B2(w1, w2, kappa, rho, eta, tau)
beta = kappa - 1i * rho * eta * w1;
D = sqrt(beta ^ 2 + eta * w1 * (1i + w1));
h = (beta - D - 1i .* w2 * eta ^ 2) / (beta + D - 1i * w2 .* eta ^ 2);

%To ease up calculations
BmD = beta - D;
hemdt = h  * exp (-D * tau);

%Return values
B2 =1 / eta * (BmD - (beta + D) * hemdt) / (1 - hemdt);
end

