function V = V_T(h, a1, b1, a2, b2, K, N)

k = 0:(N - 1);
chi = chi_coef(k, 0, h, a1, b1);
psi_1 = psi_coef(k, 0, h, a1, b1);
psi_2 = psi_coef(k, a2, b2, a2, b2);

V = (chi - psi_1)' * psi_2; % compute in vector form

%Finally
V = 2/(b1 - a1) * 2/(b2-a2) * K * V;

end