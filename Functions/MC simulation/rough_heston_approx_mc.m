function [S, V] = rough_heston_approx_mc(S0, V0, rho, kappa, theta, T, r, nu, H, npath, N)
% Uses the discritation scheme as presented in Emre Erkan 2020.
%
%  Usage:      [S, V] = heston_approx_mc(..);
%               
%  Inputs:      S0      Intial price
%               V0      Initial volatility
%               rho     correlation between S and V
%               kappa   mean-reversion speed of V
%               theta   mean of V
%               T       time to maturity
%               r       riskfree interest rate
%               nu      .. 
%               H       roughness parameter
%               npath  number of paths to be simulated
%               N      number of steps 
%  Output:      S      N x npath matrix of simulated prices
%               V      N x npath matrix of simulated volatilities

tic()

%Transformed parameters
dt = T/N;
alpha = H + 0.5; 


i_ = 0:(N - 1);
d1 = (N - i_) .^ alpha - (N - i_ - 1) .^ alpha; 
d2 = (N - i_) .^ (alpha - 1) - (N - i_ - 1) .^ (alpha - 1);
d2(N) = -1;


S = S0 * ones(N + 1, npath);
V = V0 * ones(N + 1, npath); 
 
B = zeros(N + 1, npath);
for i = 1:N 
    z = randn(2, npath / 2);
    Z = [z, -z]; %antithetic samplic
    
    %Update B
    B(i + 1, :) = B(i,:) + (rho * Z(1, :) + sqrt(1 - rho ^ 2) * Z(2, :)) * sqrt(dt);

    %Update S
    dS = r * S(i, :) * dt + sqrt(V(i, :)) .* S(i, :)* sqrt(dt) .* Z(1, :);
    S(i + 1,:) = S(i,:) + dS;
    
    %Update V
    ki = kappa * (theta - V(i, :));
    term1 = (dt ^ alpha) / gamma(alpha + 1) * d1(i) .* ki;

    li = kappa * nu * sqrt(V(i, :)) .*  B(i, :);
    term2 = (dt ^ (alpha - 1)) / (gamma(alpha) * (alpha - 1)) * d2(i) .* li;
    V(i + 1, :) = max(V(i, :) + term1 + term2, 0);
end

toc()