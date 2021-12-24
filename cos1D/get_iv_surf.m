function [iv] = get_iv_surf(model, S0, v0, K, r, q, T, eta, theta, rho, kappa, N, H)
if (model == "h"), H = 0.5;  end

[price, iv] = deal(NaN(length(T), length(K)));
for t = 1:length(T)
    [price(t, :), iv(t, :)] = cos_1d(model, S0, v0, K, r, q, T(t), eta, theta, rho, kappa, N, H);
end
end