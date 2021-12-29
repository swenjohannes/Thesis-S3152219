function phi = cfRoughHeston(w,r,v0,kappa,theta,eta,rho,H,T,nstep,X)
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

f1 = -0.5*(w.^2 + 1i*w);
func = @(x) f1 + (1i*rho*eta * w - kappa) .* x + 0.5 * eta^2 * x.^2;

if (~ exist ( 'X' , 'var' )), X = 0; end 
if (~ exist ( 'nstep' , 'var' )), nstep = 2*length(w); end 

tgrid = linspace(0,T,nstep+1);
dt = tgrid(2)-tgrid(1);
alpha = H + 0.5;
DaGa2 = dt^alpha / gamma(alpha +2);
%k = (0:nstep); 
h = zeros(length(w),nstep+1);
Dh = h;
%A = zeros(nstep+1);
%A(1,:) = DaGa2 * ((k-1).^(alpha+1) - (k-1-alpha).*k.^alpha);

Dh(:,1) = func(h(:,1));

ak = dt^alpha / gamma(alpha + 2);
j = (1:nstep);
b = dt^alpha / gamma(alpha + 1) * ((nstep+1 - j).^alpha - (nstep+1 - j - 1).^alpha);
aa = dt^alpha / gamma(alpha +2) * ((nstep+1 - j + 1).^(alpha+1) + (nstep+1 - j - 1).^(alpha+1) - 2*(nstep+1-j).^(alpha+1)); 
for k = 2:nstep+1
    %j = (1:k-1);
    %b = dt^alpha / gamma(alpha + 1) * ((k - j).^alpha - (k - j - 1).^alpha); 
    hp = Dh(:,1:k-1) * b(end-k+2:end)';  
%    figure(3)
%    plot(abs(hp))
%    title('hp')

    %a = dt^alpha / gamma(alpha +2) * ((k - j + 1).^(alpha+1) + (k - j - 1).^(alpha+1) - 2*(k-j).^(alpha+1)); 
    a = aa(end-k+2:end);
    a(1) = ak * ((k-2)^(alpha+1) - (k-2-alpha)*(k-1)^alpha);
    %a(k) = dt^alpha / gamma(alpha + 2);

    h(:,k) = Dh(:,1:k-1) * a' + ak * func(hp); 
    Dh(:,k) = func(h(:,k));
%     figure(1)
%     surf(abs(h))
%     title('h')
%     figure(2)
%     surf(abs(Dh))
%     title('Dh')
%     zyx = 1;
end

phi = exp( 1i * w * r * T ...
          + theta*kappa*sum((h(:, 1:end-1) + h(:,2:end))/2, 2) *dt ...
          + v0*sum((Dh(:, 1:end-1) + Dh(:,2:end))/2, 2) *dt ...
          + 1i * w * X);

%     function y = func(x)
%         y1 = 0.5*(w.^2 + 1i*w);
%         y2 = (1i*rho*eta * w - kappa) .* x;
%         y3 = 0.5 * eta^2 * x.^2;
%         figure(4)
%         plot(w,abs(x),w,abs(y1),w,abs(y2),w,abs(y3));
%         legend('x','y1','y2','y3');
%         y = -y1 + y2 + y3;
%     end
end

% c1 = - 0.5 *(w .^2 + 1i * w);
% beta = kappa - 1i * rho * eta * w; 
% c3 = 0.5*eta ^ 2;
% f = @(y) (c1 - beta .*y + c3*y.^2);
% 
% [psi,Dalpha_psi] = SolveVIE(f, 0, alpha,T, nstep, length(w)); %solve VIE SolveVIE(f, y0, alpha, T, N, M)
% 
% phi2 = exp( 1i * w * r * T ...
%           + theta*kappa*sum((psi(:, 1:end-1) + psi(:,2:end))/2, 2) *dt...
%            + v0.*sum((Dalpha_psi(:, 1:end-1) + Dalpha_psi(:,2:end))/2, 2) *dt);
% 
% 
% max(abs(psi-h)./abs(h+1e-14))
% 
% max(abs(Dh-Dalpha_psi)./abs(Dh))
% 
% max(abs(phi-phi2))
% 

% function f = func(u,x)
% f = -0.5*(u.^2 + 1i*u) + (1i*rho*eta * u - kappa) .* x + 0.5 * eta^2 * x.^2;
% end
% end
