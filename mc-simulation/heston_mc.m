function [S, V] = heston_mc(S0, v0, rho, kappa, theta, T, r, q, eta, npath, N)
%{
 Description: MC simulation of the Heston model using an euler scheme
              
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
 References:

%}

dt = T/N;                               %Time delta
S = S0 * ones(N + 1, npath);            %Matrix to store stock results
V = v0 * ones(N + 1, npath);            %Matrix to store variance results

for t = 1:N
    z =  randn(1, npath/2);                             %antithetic sampling
    Zs = [z, -z];                                       %stock errors
    Zv = rho * Zs + sqrt(1 - rho ^2) * randn(1, npath); %variance error

    S(t + 1, :) = S(t, :) + (r - q) * S(t, :) * dt + ... 
                                sqrt(V(t, :) * dt) .* S(t, :) .* Zs;
    

    V(t + 1, : ) = V(t, :) + kappa * (theta - V(t, :)) * dt ...
                            + eta * sqrt(V(t, :) * dt) .* Zv;
    V(t + 1, :) = max(V(t + 1, :), 0);    %variance can't become negative
end
