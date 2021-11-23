function V = V_T(h, a, b, K, N)

k = 0:(N - 1);
chi = chi_coef(k, 0, h, a, b);
psi_1 = psi_coef(k, 0, h, a, b);
psi_2 = psi_coef(k, a, b, a, b);

V = (chi - psi_1)' * psi_2; % compute in vector form

%Finally
V = (2 / (b - a)) ^ 2 * K * V;

end