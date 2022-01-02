function [phi,B2] = cfRoughHeston2(w,r,v0,kappa,theta,eta,rho,H,T,nstep,X)
% clear all
% 
% w = (0:12)'/12*pi;
% 
% v0 = 0.2^2;
% kappa = 1.7;
% theta = 0.22;
% eta = 0.5;
% rho = -0.7;  
% H = 0.1;
% T = 1;
% nstep = 10;
% r = 0.05;

if (~ exist ( 'X' , 'var' )), X = 0; end 
if (~ exist ( 'nstep' , 'var' )), nstep = 30; end % TODO: the optimal nstep is still to be found

alpha = H + 0.5;
f1 = 0.5*(w.^2 + 1i*w);
func = @(x) -f1 + (1i*rho*eta * w - kappa) .* x + 0.5 * eta^2 * x.^2;

tgrid = linspace(0,T,nstep+1);
dt = tgrid(2)-tgrid(1);
h = zeros(length(w),nstep+1);
Dh = h;

Dh(:,1) = func(h(:,1));

ak = dt^alpha / gamma(alpha + 2);
j = (1:nstep);
aa = ak * ((nstep+1 - j + 1).^(alpha+1) + (nstep+1 - j - 1).^(alpha+1) - 2*(nstep+1-j).^(alpha+1)); 

for k = 2:nstep+1
    a = aa(end-k+2:end);
    a(1) = ak * ((k-2)^(alpha+1) - (k-2-alpha)*(k-1)^alpha);
    
    A = ak * 0.5 * eta^2;
    B = ak * (1i*rho*eta * w - kappa) - 1;
    C = Dh(:,1:k-1) * a' - ak * f1;
    D = (B.^2 - 4*A.*C) ./ A.^2 /4;
    
    nn = 1;
    h(:,k) = -B./A/2 + sqrt(abs(D)).*exp(1i*angle(D)/2 + 1i * pi * nn); 
    Dh(:,k) = func(h(:,k));
end

B2 = sum((Dh(:, 1:end-1) + Dh(:,2:end))/2, 2) *dt;
phi = exp( 1i * w * r * T ...
          + theta*kappa*sum((h(:, 1:end-1) + h(:,2:end))/2, 2) *dt ...
          + v0*sum((Dh(:, 1:end-1) + Dh(:,2:end))/2, 2) *dt ...
          + 1i * w * X);
end

