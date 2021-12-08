function [S, V] = rough_heston_euler_scheme(S0, V0, rho, kappa, theta, T, nu, H, npath, N)
% Uses the discritation scheme as presented in Emre Erkan 2020.
%
%  Usage:      [S, V] = rough_heston_euler_scheme(..);
%               
%  Inputs:      S0      Intial price
%               V0      Initial volatility
%               rho     correlation between S and V
%               kappa   mean-reversion speed of V
%               theta   mean of V
%               T       time to maturity
%               nu      .. 
%               H       roughness parameter
%               npath  number of paths to be simulated
%               N      number of steps 
%  Output:      S      N x npath matrix of simulated prices
%               V      N x npath matrix of simulated volatilities

tic()

%Transformed parameters
dt = T/N;
C = [1,rho; rho, 1];
L = chol(C);

Y = log(S0) * ones(N + 1, npath); %Matrix of simulated prices
V = V0 * ones(N + 1, npath);      %Matrix of simulated vol.

dY = zeros(N, npath); %To keep track of dy
W = zeros(N, npath);  %To keep track of W

K = @(t) 1 / gamma(H) * t ^ (H - 1 / 2);

for k = 2:(N + 1)
    z = randn(2, npath / 2);
    z = [z, -z]; %antithetic samplic
    Z = L * z;  

    %Update Y
    dY(k - 1, :) = - 1 / 2 * V(k - 1, :) * dt + sqrt(V(k - 1, :)) .* Z(1, :) * sqrt(dt);
    Y(k, :) = Y(k - 1, :) + dY(k - 1, :);

    %Update V
    W(k - 1, :) = Z(2, :) * sqrt(dt); %store the variance error process!
    
    dV = 0;
    for m = 1:(k-1)
        dV = dV + K((k-m) * dt) * kappa * (theta - V(m, :)) * dt + ...
        K((k-m) * dt) * nu * sqrt(V(m, :)) .* W(m, :);
    end
    V(k, :) = max(V(1, :) + dV, 0);
    %% 
 
end

S = exp(Y);

toc()