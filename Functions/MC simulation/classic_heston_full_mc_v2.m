function [S, V] = classic_heston_full_mc(S0, V0, rho, kappa, theta, T, r, nu, npath, N)
% Similar to the lecture slides (quantitative finance)
%
%  Usage:      [S, V] = classic_heston_full_mc(..);
%               
%  Inputs:      S0      Intial price
%               V0      Initial volatility
%               rho     correlation between S and V
%               kappa   mean-reversion speed of V
%               theta   mean of V
%               T       time to maturity
%               r       riskfree interest rate
%               nu      .. 
%               npath  number of paths to be simulated
%               N      number of steps 
%  Output:      S      N x npath matrix of simulated prices
%               V      N x npath matrix of simulated volatilities

tic()

%Transformed parameters
dt = T/N;


%Simulate random draws:
S = S0 * ones(N + 1, npath);
V = V0 * ones(N + 1, npath); 

for t = 1:N
    Zs = randn(1, npath);
    Zv = rho * Zs + sqrt(1 - rho ^2) * randn(1, npath);

    S(t + 1, :) = S(t, :) .* exp(sqrt(V(t, :) * dt) .* Zs ... 
                            - V(t, :) * dt  / 2);
    V(t + 1, : ) = V(t, :) + kappa * (theta - V(t, :)) * dt ...
        + nu * sqrt(V(t, :) * dt) .* Zv;

    V(t + 1, :) = max(V(t + 1, :), 0);
end

toc()