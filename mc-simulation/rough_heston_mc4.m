function [S, V, time] = rough_heston_mc4(S0, v0, rho, kappa, theta, T, r, q, eta, H, npath, N)
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

% N = 10, npath = 6

tic();
alpha = H + 1/2;
dt = T/N;                       %Time steps
S = NaN(N + 1, npath);          %Matrix to store stock results
V = NaN(N + 1, npath);          %Matrix to store variance results

S(1, :) = S0;
V(1, :) = v0;

W = NaN(N, npath);                        %To keep track of W

t_obs = (1:N) * dt;                             %make column vector
K = 1 / gamma(alpha) * t_obs .^ (alpha - 1);    %Calculate kernel values
K = flip(K);

for t = 1:N
    z =  randn(1, npath/2);                             %antithetic sampling
    Zs = [z, -z];                                       %stock errors
    Zv = rho * Zs + sqrt(1 - rho ^2) * randn(1, npath); %variance error

    S(t + 1, :) = S(t, :) + (r - q) * S(t, :) * dt + ... 
                                sqrt(V(t, :) * dt) .* S(t, :) .* Zs;

    W(t, :) = sqrt(V(t, :) * dt) .* Zv; %store the variance error process!
    
    V(t + 1, : ) = max(v0 ...
                        + K(end - t +1:end) * kappa * (theta - V(1:t, :)) * dt ...
                        + K(end - t +1:end) * eta * W(1:t, :), 0);
end

fprintf('Rough Heston MC: ')   
time = toc()