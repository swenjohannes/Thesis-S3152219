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

debug = false;

if (~ exist ( 'X' , 'var' )), X = 0; end 
if (~ exist ( 'nstep' , 'var' )), nstep = 100; end 

% if T < 14/365
%     nstep = 256;
% elseif T < 0.5
%     nstep = 512;
% else
%     nstep = 5000;
% end

alpha = H + 0.5;
f1 = 0.5*(w.^2 + 1i*w);
func = @(x) -f1 + (1i*rho*eta * w - kappa) .* x + 0.5 * eta^2 * x.^2;

% tgrid = linspace(0,T,nstep+1);
% dt = tgrid(2)-tgrid(1);
% h = zeros(length(w),nstep+1);
% Dh = h;
% 
% 
% Dh(:,1) = func(h(:,1));
% 
% ak = dt^alpha / gamma(alpha + 2);
% j = (1:nstep);
% b = dt^alpha / gamma(alpha + 1) * ((nstep+1 - j).^alpha - (nstep+1 - j - 1).^alpha);
% aa = dt^alpha / gamma(alpha +2) * ((nstep+1 - j + 1).^(alpha+1) + (nstep+1 - j - 1).^(alpha+1) - 2*(nstep+1-j).^(alpha+1)); 
% for k = 2:nstep+1
%     %j = (1:k-1);
%     %b = dt^alpha / gamma(alpha + 1) * ((k - j).^alpha - (k - j - 1).^alpha); 
% %     hp = Dh(:,1:k-1) * b(end-k+2:end)';  
% 
%     %a = dt^alpha / gamma(alpha +2) * ((k - j + 1).^(alpha+1) + (k - j - 1).^(alpha+1) - 2*(k-j).^(alpha+1)); 
%     a = aa(end-k+2:end);
%     a(1) = ak * ((k-2)^(alpha+1) - (k-2-alpha)*(k-1)^alpha);
%     %a(k) = dt^alpha / gamma(alpha + 2);
%     hp = Dh(:,1:k-1) * a' + ak * func(h(:,k-1));
%     %hp = Dh(:,1:k-1) * b(end-k+2:end)';
% 
%     if debug    
%         figure(3)
%         plot(abs(hp))
%         title('hp')
%     end
%     h(:,k) = Dh(:,1:k-1) * a' + ak * func(hp); 
%     Dh(:,k) = func(h(:,k));
%     
%     if debug
%         figure(1)
%         plot(1:k,abs(h(end-2:end,1:k)))
%         %view([180 0])
%         title('h')
%         figure(2)
%         surf(abs(Dh))
%         view([180 0])
%         title('Dh')
%         zyx = 1;
%     end
% end
h = nan;
while any(isnan(h(:,end))) && (nstep < 20000)
    [h,Dh] = volterra(nstep);
    nstep = nstep*2;
end

dt = T/(size(h,2)-1);
phi = exp( 1i * w * r * T ...
          + theta*kappa*sum((h(:, 1:end-1) + h(:,2:end))/2, 2) *dt ...
          + v0*sum((Dh(:, 1:end-1) + Dh(:,2:end))/2, 2) *dt ...
          + 1i * w * X);

    function [g,Dg] = volterra(nstepin)
        tgrid = linspace(0,T,nstepin+1);
        dt = tgrid(2)-tgrid(1);
        g = zeros(length(w),nstepin+1);
        Dg = g;
        
        Dg(:,1) = func(g(:,1));
        ak = dt^alpha / gamma(alpha + 2);
        j = (1:nstepin);
        aa = dt^alpha / gamma(alpha +2) * ((nstepin+1 - j + 1).^(alpha+1) + (nstepin+1 - j - 1).^(alpha+1) - 2*(nstepin+1-j).^(alpha+1)); 

        for k = 2:nstepin+1
            a = aa(end-k+2:end);
            a(1) = ak * ((k-2)^(alpha+1) - (k-2-alpha)*(k-1)^alpha);
            gp = Dg(:,1:k-1) * a' + ak * func(g(:,k-1));
            g(:,k) = Dg(:,1:k-1) * a' + ak * func(gp); 
            Dg(:,k) = func(g(:,k));
            if any(isinf(Dg(:,k)))
                g = nan(length(w),nstepin+1);
                %warning(['Volterra intergarl diverges at step ' int2str(k) '. The number of steps is doubled.'])
                break
            end
        end        
    end
%     function y = func(x)
%         y1 = 0.5*(w.^2 + 1i*w);
%         y2 = (1i*rho*eta * w - kappa) .* x;
%         y3 = 0.5 * eta^2 * x.^2;
%         if debug
%             figure(4)
%             plot(w,abs(x),w,abs(y1),w,abs(y2),w,abs(y3));
%             legend('x','y1','y2','y3');
%             zyx = 1;
%         end
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
