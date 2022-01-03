function [phiA, B2] = phi_A_B2_rough_heston2(u1,u2, r,kappa,theta,eta,rho,H,T,nstep)

%u1 = (0:50)'*pi / (b1 - a1);
%u2 = (0:50)'*pi / (b2 - a2);
% v0 = 0.2^2;
% kappa = 1.7;
% theta = 0.22;
% eta = 0.5;
% rho = -0.7;  
% H = 0.1;
% T = 0.5;
% nstep = 100;
% r = 0.05;

if (~ exist ( 'nstep' , 'var' )), nstep = 30; end % TODO: the optimal nstep is still to be found

alpha = H + 0.5;
[phiA,B2] = deal(zeros(length(u1),length(u2)));

f1 = 0.5*(u1.^2 + 1i*u1);
func = @(x) -f1 + (1i*rho*eta * u1 - kappa) .* x + 0.5 * eta^2 * x.^2;

tgrid = linspace(0,T,nstep+1);
dt = tgrid(2)-tgrid(1);

for m = 1:length(u2)
    %m = 5;
    %h = zeros(length(u1),nstep+1);
    h = 1i*u2(m) * ones(length(u1),nstep+1);
    Dh = zeros(length(u1),nstep+1);

    Dh(:,1) = func(h(:,1));

    ak = dt^alpha / gamma(alpha + 2);
    j = (1:nstep);
    aa = ak * ((nstep+1 - j + 1).^(alpha+1) + (nstep+1 - j - 1).^(alpha+1) - 2*(nstep+1-j).^(alpha+1)); 

    for k = 2:nstep+1
        a = aa(end-k+2:end);
        a(1) = ak * ((k-2)^(alpha+1) - (k-2-alpha)*(k-1)^alpha);

        A = ak * 0.5 * eta^2;
        B = ak * (1i*rho*eta * u1 - kappa) - 1;
        C = 1i * u2(m) * tgrid(k)^(alpha-1)/gamma(alpha) + Dh(:,1:k-1) * a' - ak * f1;
        D = (B.^2 - 4*A.*C) ./ A.^2 /4;

        nn = 1;
        h(:,k) = -B./A/2 + sqrt(abs(D)).*exp(1i*angle(D)/2 + 1i * pi * nn); 
        Dh(:,k) = func(h(:,k));
    end
    %plot_imag(h)
    %plot_imag(Dh)
    phiA(:,m) = exp(1i * r * T * u1 + theta*kappa * sum((h(:, 1:end-1) + h(:,2:end))/2, 2) *dt);
    B2(:,m) = sum((Dh(:, 1:end-1) + Dh(:,2:end))/2, 2) *dt + 1i*u2(m);
end

