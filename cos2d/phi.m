function phi = phi(x, v, w1, w2, kappa, rho, eta, theta, r, tau)

%Get A, B2
ePhiA = phi_A(w1, w2, kappa, rho, eta, theta, r, tau);
B = B2(w1, w2, kappa, rho, eta, tau);
phi = exp(1i * w1 * x) * exp(B * v) * ePhiA;
end