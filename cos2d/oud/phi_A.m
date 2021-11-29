function phi = phi_A(w1, w2, kappa, rho, eta, theta, r, tau)

beta = kappa - 1i * rho * eta * w1;
D = sqrt(beta ^ 2 + eta * w1 * (1i + w1));
h = (beta - D - 1i * w2 * eta ^ 2) / (beta + D - 1i * w2 * eta ^ 2);

%To ease up calculations
BmD = beta - D;
H = h  * exp (-D * tau);

%Return values
B2 =1 / eta * (BmD - (beta + D) * H) / (1 - H);
A  = 1i * w1 * r * tau + ...
    kappa * theta / eta ^ 2 * (BmD * tau - 2 * log((H - 1) / (h-1)));
phi = exp(A); 

end