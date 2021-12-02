function V = mc_V(V0, kappa, theta, T, nu, npath, N)
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

%Transformed parameters
dt = T/N;

%Simulate random draws:
V = V0 * ones([1, npath]);

for t = 1:N
    Zv = randn(1, npath);
    V = V + kappa * (theta - V) * dt ...
        + nu * sqrt(V * dt) .* Zv;
    V = max(V, 0);
end

