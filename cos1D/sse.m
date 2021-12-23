function [val] = sse(vol_dat, S0, v0, K, r, q, T, eta, theta, rho, kappa)
val = sum((vol_dat - get_iv_surf(S0, v0, K, r, q, T, eta, theta, rho, kappa)).^2, 'all');
end