function y = cf2dRoughHeston_mc(u1,u2, x0,v0, r,kappa,theta,eta,rho,H,T, nstep,npath)

if (~ exist ( 'nstep' , 'var' )), nstep = 1e3; end
if (~ exist ( 'npath' , 'var' )), npath = 1e3; end

u1 = u1(:);

tgrid = linspace(0,T,nstep+1);
dt = tgrid(2) - tgrid(1);               %Time delta
alpha = H+0.5;
Kt = 1/gamma(alpha) * (tgrid(2:end).^(alpha-1))';
M = zeros(npath, nstep+1);              %Increments of the variance

Y = x0*ones(npath, nstep+1);            %Log price
V = v0*ones(npath, nstep+1);            %Matrix to store variance results
C = [1 rho; rho 1];
U = chol(C);

for i = 1:nstep
  Z1 = randn(npath/2,2);
  Z = [Z1;-Z1];
  Z = Z * U;
  dY = (r - 0.5*V(:,i))*dt + sqrt(V(:,i)) .*Z(:,1)*sqrt(dt);
  Y(:,i+1) = Y(:,i) + dY;
    
  M(:,i) = sqrt(V(:,i)).*Z(:,2) * (sqrt(dt)*eta);
  Kts = flipud(Kt(1:i));
  dV = ((theta - V(:,1:i)) *(kappa*dt) + M(:,1:i)) * Kts;
  V(:,i+1) = max(v0 + dV,0);
end
XT = Y(:,end);
VT = V(:,end);

y = zeros(length(u1),length(u2));

%A = exp((1i * u1) * XT');
for m = 1:length(u2)
    %y(:,m) = mean(A .* exp((1i * u2(m)) * ones(size(u1)) * VT'),2);
    y(:,m) = mean(exp((1i * u1) * XT' + (1i * u2(m)) * ones(size(u1)) * VT'),2);    
end
end