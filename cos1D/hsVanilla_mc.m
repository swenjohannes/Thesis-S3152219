function price = hsVanilla_mc(K,pc,T, S0,r,q,v0,kappa,theta,eta,rho, npath,nstep)

dt = T/nstep;                               %Time delta
Y = log(S0)*ones(npath, nstep+1);       %Log price
V = v0*ones(npath, nstep+1);            %Matrix to store variance results
C = [1 rho; rho 1];
U = chol(C);

for i = 1:nstep
  Z1 = randn(npath/2,2);
  Z = [Z1;-Z1];
  Z = Z * U;
  dY = (r-q - 0.5*V(:,i))*dt + sqrt(V(:,i)) .*Z(:,1)*sqrt(dt);
  Y(:,i+1) = Y(:,i) + dY;
  
  dV = (theta - V(:,i))*(kappa*dt) + sqrt(V(:,i)).*Z(:,2)*(eta*sqrt(dt));
  V(:,i+1) = max(V(:,i) + dV,0);  
end
%S = exp(Y);
ST = exp(Y(:,end));
%Fmc = mean(ST)

price = zeros(size(K));
if pc >= 0 % Call
    for i = 1:length(K)
      price(i) = exp(-r*T)*mean((ST-K(i)).*(ST>K(i)));
    end
else % Put
    for i = 1:length(K)
      price(i) = exp(-r*T)*mean((K(i)-ST).*(ST<K(i)));
    end
end
end