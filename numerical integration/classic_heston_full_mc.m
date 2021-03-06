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
C = [1,rho; rho, 1];
L = chol(C);

%Simulate random draws:
S = S0 * ones(N + 1, npath);
V = V0 * ones(N + 1, npath); 

for i = 1:N
    z = randn(2, npath / 2);
    z = [z, -z]; %antithetic samplic
    Z = L * z;  

    dS = r * S(i, :) * dt + sqrt(V(i, :)) .* S(i, :)* sqrt(dt) .* Z(1, :);
    S(i+1,:) = S(i,:) + dS;
    dV = kappa* (theta - V(i, :)) * dt + nu * sqrt(V(i, :)) * sqrt(dt) .* Z(2, :);
    V(i+1, :) = max(V(i,:) + dV, 0);
end

toc()