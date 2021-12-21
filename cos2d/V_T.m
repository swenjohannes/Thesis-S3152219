function V = V_T(h, a1, b1, a2, b2, K, N)

k = 0:(N - 1);
chi_k1 = chi_coef(k, 0, h, a1, b1);
psi_k1 = psi_coef(k, 0, h, a1, b1);
psi_k2 = psi_coef(k, a2, b2, a2, b2);

V = (chi_k1 - psi_k1)' * psi_k2; % compute in vector form

%Finally
omega1 = 2 / (b1 - a1);
omega2 = 2 / (b2 - a2);
V = omega1 * omega2 * K * V;
end