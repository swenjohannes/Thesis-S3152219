function y = cf2dHeston_mc(u1,u2, x0,v0, r,kappa,theta,eta,rho,T,nstep,npath)

if (~ exist ( 'nstep' , 'var' )), nstep = 1e5; end
if (~ exist ( 'npath' , 'var' )), npath = 1e3; end

u1 = u1(:);

dt = T/nstep;                               %Time delta
Y = x0*ones(npath, nstep+1);       %Log price
V = v0*ones(npath, nstep+1);            %Matrix to store variance results
C = [1 rho; rho 1];
U = chol(C);

for i = 1:nstep
  Z1 = randn(npath/2,2);
  Z = [Z1;-Z1];
  Z = Z * U;
  dY = (r - 0.5*V(:,i))*dt + sqrt(V(:,i)) .*Z(:,1)*sqrt(dt);
  Y(:,i+1) = Y(:,i) + dY;

  dV = (theta - V(:,i))*(kappa*dt) + sqrt(V(:,i)).*Z(:,2)*(eta*sqrt(dt));
  V(:,i+1) = max(V(:,i) + dV,0);  
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