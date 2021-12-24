function [S_T, time] = heston_mc_final_only(S, V, rho, kappa, theta, T, r, q, eta, npath, N, Nobs)
%{
 Description: MC simulation of the Heston model using an euler scheme, only
              returns final value!
              
 Inputs:      S      Intial price
              V      Initial volatility              

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
tic();
dt = T/N;                               %Time delta
t_obs = ((1:Nobs) * N / Nobs) + 1;  %index of observation moment
S_T = NaN([length(t_obs), npath]);

for t = 2:(N +1)
    z =  randn(1, npath/2);                             %antithetic sampling
    Zs = [z, -z];                                       %stock errors
    Zv = rho * Zs + sqrt(1 - rho ^2) * randn(1, npath); %variance error

    S = S + (r - q) * S * dt + sqrt(V * dt) .* S .* Zs;
    if sum(t == t_obs) == 1
        S_T(find(t == t_obs ), :) = S; %check if we need to store te result
    end
    V = V  + kappa * (theta - V) * dt ...
                            + eta * sqrt(V * dt) .* Zv;
    V = max(V, 0);                  %variance can't become negative
end
fprintf('Heston MC: ')   
time = toc()