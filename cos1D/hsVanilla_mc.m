function price = hsVanilla_mc(K,pc,T, S0,r,q,v0,kappa,theta,eta,rho, npath,nstep)

dt = T/nstep;                               %Time delta
Y = log(S0)*ones(npath, nstep+1);       %Log price
V = v0*ones(npath, nstep+1);            %Matrix to store variance results
Ybs = log(S0)*ones(npath, nstep+1);     %Log price in BS model
Vbs = sqrt(v0);                         %Volatility in BS model
C = [1 rho; rho 1];
U = chol(C);

for i = 1:nstep
  Z1 = randn(npath/2,2);
  Z = [Z1;-Z1];
  Z = Z * U;
  dY = (r-q - 0.5*V(:,i))*dt + sqrt(V(:,i)) .*Z(:,1)*sqrt(dt);
  Y(:,i+1) = Y(:,i) + dY;

  dYbs = (r-q - 0.5*Vbs^2)*dt + Vbs .*Z(:,1)*sqrt(dt);
  Ybs(:,i+1) = Ybs(:,i) + dYbs;
  
  dV = (theta - V(:,i))*(kappa*dt) + sqrt(V(:,i)).*Z(:,2)*(eta*sqrt(dt));
  V(:,i+1) = max(V(:,i) + dV,0);  
end
%S = exp(Y);
ST = exp(Y(:,end));
STbs = exp(Ybs(:,end));
%Fmc = mean(ST)

price = zeros(size(K));
if pc >= 0 % Call
    for i = 1:length(K)
      price(i) = exp(-r*T)*mean((ST-K(i)).*(ST>K(i))) - ( exp(-r*T)*mean((STbs-K(i)).*(STbs>K(i))) - bsprice(S0,K(i),r,q,Vbs,T,pc));
    end
else % Put
    for i = 1:length(K)
      price(i) = exp(-r*T)*mean((K(i)-ST).*(ST<K(i)))- ( exp(-r*T)*mean((K(i)-STbs).*(STbs < K(i))) - bsprice(S0,K(i),r,q,Vbs,T,pc));
    end
end
end