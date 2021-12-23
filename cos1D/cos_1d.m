function [price, iv] = cos_1d(model, S0, v0, K, r, q, T, eta, theta, rho, kappa, N, H)
x = log(S0 ./ K);

k = (0:(N -1))';

if model == "h"
    %Obtain trucation ranges from cumulants
    [c1, c2, ~]  = heston_cumulants_v1(r, q, kappa, theta, v0, eta, rho, T);
    [a, b] = cos_truncation_range_v2(c1,c2,0,12); %Obtain a and b from cumulants
    phi = chfun_heston(v0, kappa, rho, eta, theta, r, q, T, k * pi / (b -a));
else
    b = 12 * sqrt(T); %Obtain a and b from L rule!
    a = -b;
    phi = chfun_rough_heston(r, q, kappa, theta, v0, eta, rho, T, k * pi / (b -a), H, 250);
end
Uk = 2 / (b - a) *(chi_coef(k, 0, b, a, b) - psi_coef(k, 0, b, a, b));

F = phi .* Uk .* exp(1i * k * pi * (x - a) / (b - a));
F(1, :) = 0.5 * F(1, :); %weight first term by a half

price = K * exp(-r * T) .* real(sum(F, 1));
if any(price < 0)
    iv = 2; %Set artifically high
else
    iv = blsimpv(S0,K,r,T,price, 'Yield',q);
end

end