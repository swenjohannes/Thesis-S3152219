function price = hrBarrier_mc(K,B,pc,du,io,T,nobs, S0,r,q,v0,kappa,theta,eta,rho,H, npath,nstep,seed)
%{
 Description: MC simulation of the Heston model using an euler scheme
              
 Inputs:
            K       strike
            B       barrier
            pc      put/call -1/+1
            du      down/up -1/+1
            io      in/out -1/+1
            T       maturity
            nobs    number of equidistant observations; 0 for continuous

            S0      Intial price
            v0      Initial volatility              
            r       riskfree interest rate
            q       divident yield
            kappa   mean-reversion speed of V
            theta   mean of V
            eta     volatility of volatility
            rho     correlation between S and V
            H       roughness parameter

            npath   number of paths to be simulated
            nstep       number of steps 
 References:

%}
if (exist ( 'seed' , 'var' )), rng(seed); end 

tgrid = linspace(0,T,nstep+1);
dt = tgrid(2) - tgrid(1);               %Time delta
alpha = H+0.5;
Kt = 1/gamma(alpha) * (tgrid(2:end).^(alpha-1))';
M = zeros(npath, nstep+1);             %Increments of the variance

Y = log(S0)*ones(npath, nstep+1);       %Log price
V = v0*ones(npath, nstep+1);            %Matrix to store variance results
I = true(npath,1);                      %Matrix for knock indicator
C = [1 rho; rho 1];
U = chol(C);
Nobs = round((1:nobs)/nobs*nstep);
logB = log(B);

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
  
  if any(Nobs == i+1)
    I = I & (Y(:,i+1) < logB);
  end
end
S = exp(Y);
ST = S(:,end);
Fmc = mean(ST);
disp(['Fmc = ' num2str(Fmc)])

price = zeros(size(K));
for i = 1:length(K)
  price(i) = exp(-r*T)*mean((ST-K(i)).*(ST>K(i)).*I);
end
