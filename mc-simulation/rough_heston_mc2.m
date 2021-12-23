function [S, V] = rough_heston_mc2(S0, v0, rho, kappa, theta, T, r, q, eta, H, npath, N)
%{
 Description: MC simulation of rough volatility using an euler scheme,
              based on ..
              
 Inputs:      S0      Intial price
              v0      Initial volatility              

              rho     correlation between S and V
              kappa   mean-reversion speed of V
              theta   mean of V
              T       time to maturity
              r       riskfree interest rate
              q       divident yield
              eta     volatility of volatility
              H       roughness parameter
              npath   number of paths to be simulated
              N       number of steps 
 Output:      S       N x npath matrix of simulated prices
              V       N x npath matrix of simulated volatilities
%}
tic()
dt = T/N;                       %Time steps
C = [1,rho; rho, 1];
L = chol(C);

Y = log(S0) * ones(N + 1, npath);           %Matrix of simulated prices
V = v0 * ones(N + 1, npath);                %Matrix of simulated vol.
dY = zeros(N, npath);                       %To keep track of dy
W = zeros(N, npath);                        %To keep track of W

t_obs = (1:N) * dt;                         %make column vector
K = 1 / gamma(H) * t_obs .^ (H - 1 / 2);    %Calculate kernel values

for k = 2:(N + 1)
    z = randn(2, npath / 2);
    z = [z, -z];                            %antithetic samplic
    Z = L * z;                              %correlation coefficient

    m = k - 1; 
    dY = (r - q - 0.5 * V(k - 1, :)) * dt + sqrt(V(k - 1, :)) .* Z(1, :) * sqrt(dt);
    Y(k, :) = Y(k - 1, :) + dY;    %Update Y
    
    W(m, :) = sqrt(V(m, :) * dt) .* Z(2, :); %store the variance error process!
    
    Km = flip(K(1:m));
    dV =  Km * (theta - V(1:m, :))  * kappa * dt  ...
             + eta * Km * W(1:m, :); %Update V 
    V(k, :) = max(V(m, :) + dV, 0);
end

S = exp(Y);
toc()