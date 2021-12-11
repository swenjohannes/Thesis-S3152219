function y = cf2dHeston(u1,u2, x,v, r,tau,kappa,nu,eta,rho)
% 2D characteristic function of the Heston model

u1 = u1(:);
u2 = u2(:);
x = x(:)';
v = v(:)';

beta = kappa - 1i * rho * u1;
D = sqrt(beta^2 + eta^2 * u1 * (1i + u1));
h = (beta - D - 1i * u2 * eta^2)./(beta + D - 1i * u2 * eta^2);
eDt = exp(-D * tau);

B1 = 1i * u1;
B2 = 1/eta^2 * (beta - D - (beta + D) .* h .* eDt) ./ (1 - h .* eDt);
A = 1i * u1 * r * tau + kappa * nu / eta^2 *((beta - D) * tau - 2 * ln((h.*eDt - 1)./(h-1)));

y = exp(B1 * x + B2 * v + A);