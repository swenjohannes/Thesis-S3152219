function psi = psi_N(v, a,  T, rho, sigma, kappa, theta, x0, v0, rd, rf)
phi_N_ = @(u) phi_N(u, T, rho, sigma, kappa, theta, x0, v0, rd, rf);

v_tilde = v - 1i * a;
n = length(v);
u = v_tilde - [repelem(0, n-1) 1i];
denom = (1i * v_tilde(n) + 1) * 1i * prod(v_tilde);

psi = exp(-rd * T) * phi_N_(u) / denom;

end
