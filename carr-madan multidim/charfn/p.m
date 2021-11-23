function p = p(u, kappa, rho, sigma)
p = (0.5 - kappa * rho / sigma - 1/2 * 1i * u * (1 - rho ^2 )) * 1i * u;
end