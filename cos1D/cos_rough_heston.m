function [price, iv] = cos_rough_heston(S0, v0, K, r, q, T, eta, theta, rho, kappa, N, H)

x = log(S0 ./ K);
k = (0:(N -1))';

b = 12 * sqrt(T);
a = -b;
phi = chfun_rough_heston(r, q, kappa, theta, v0, eta, rho, T, k * pi / (b -a), H, 200);

Uk = 2 / (b - a) *(chi_coef(k, 0, b, a, b) - psi_coef(k, 0, b, a, b));

F = phi .* Uk .* exp(1i * k * pi * (x - a) / (b - a));
F(1, :) = 0.5 * F(1, :); %weight first term by a half

price = K * exp(-r * T) .* real(sum(F, 1));

price(price<0) = 0;
iv = blsimpv(S0,K,r,T,price, 'Yield',q);
end

