function price = hrVanilla_mc(K,pc,T, S0,r,q,v0,kappa,theta,eta,rho,H, npath,nstep)

tgrid = linspace(0,T,nstep+1);
dt = tgrid(2) - tgrid(1);               %Time delta
alpha = H+0.5;
Kt = 1/gamma(alpha) * (tgrid(2:end).^(alpha-1))';
M = zeros(npath, nstep+1);             %Increments of the variance

Y = log(S0)*ones(npath, nstep+1);       %Log price
Yhs = log(S0)*ones(npath, nstep+1);     %Log price in standard Heston model
Vhs = v0*ones(npath, nstep+1);          %Volatility in standard Heston model
V = v0*ones(npath, nstep+1);            %Matrix to store variance results
C = [1 rho; rho 1];
U = chol(C);

for i = 1:nstep
  Z1 = randn(npath/2,2);
  Z = [Z1;-Z1];
  Z = Z * U;
  dY = (r-q - 0.5*V(:,i))*dt + sqrt(V(:,i)) .*Z(:,1)*sqrt(dt);
  Y(:,i+1) = Y(:,i) + dY;
    
  M(:,i) = sqrt(V(:,i)).*Z(:,2)*sqrt(dt)*eta;
  Kts = flipud(Kt(1:i));
  dV = ((theta - V(:,1:i)) *(kappa*dt) + M(:,1:i)) * Kts;
  V(:,i+1) = max(v0 + dV,0);
  
  dYhs = (r-q - 0.5*Vhs(:,i))*dt + sqrt(Vhs(:,i)) .*Z(:,1)*sqrt(dt);
  Yhs(:,i+1) = Yhs(:,i) + dYhs;
  dVhs = (theta - Vhs(:,i))*(kappa*dt) + sqrt(Vhs(:,i)).*Z(:,2)*(eta*sqrt(dt));
  Vhs(:,i+1) = max(Vhs(:,i) + dVhs,0);  
end
%S = exp(Y);
ST = exp(Y(:,end));
SThs = exp(Yhs(:,end));
%Fmc = mean(ST)

price = zeros(size(K));
if pc >= 0 % Call
    for i = 1:length(K)
      price(i) = exp(-r*T)*mean((ST-K(i)).*(ST>K(i))) - ( exp(-r*T)*mean((SThs-K(i)).*(SThs>K(i))) - hsVanilla_cos(K(i),pc,T, S0,r,q,v0,kappa,theta,eta,rho));
    end
else % Put
    for i = 1:length(K)
      price(i) = exp(-r*T)*mean((K(i)-ST).*(ST<K(i))) - ( exp(-r*T)*mean((K(i)-SThs).*(SThs < K(i))) - hsVanilla_cos(K(i),pc,T, S0,r,q,v0,kappa,theta,eta,rho));
    end
end
end