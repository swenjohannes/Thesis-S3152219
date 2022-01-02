function [y,A,B1,B2] = cf2dHeston(u1,u2, x0,v0, r,kappa,theta,eta,rho,T)
% 2D characteristic function of the Heston model

u1 = u1(:);

[y,A,B1,B2] = deal(zeros(length(u1),length(u2)));

for m = 1:length(u2)
    beta = kappa - 1i * rho * eta * u1;
    d = beta.^2 + eta^2 * u1 .* (1i + u1);
    nn = 0;
    D = sqrt(abs(d)).*exp(0.5i * angle(d) + 1i*pi * nn) ;
    h = (beta - D - 1i * eta^2 * u2(m) )./(beta + D - 1i * eta^2 * u2(m));
    heDt = h.*exp(-D * T);

    B1(:,m) = 1i * u1 ;
    B2(:,m) = 1/eta^2 * (beta - D - (beta + D) .* heDt) ./ (1 - heDt);
    A(:,m) = (1i * r * T) * u1 + kappa * theta / eta^2 *((beta - D) * T - 2 * log((heDt - 1)./(h-1)));

    y = exp(B1 * x0 + B2 * v0 + A);
end