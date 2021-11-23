function [S] = bs_mc(S0, sigma, r, T, npath, N)

%Transformed parameters
dt = T/N;
S = S0 * ones(N + 1, npath);
for t = 1:N
    Zs = randn(1, npath);
    dS = S(t, :) * r * dt +  sigma  * S(t, :) .* Zs * sqrt(dt);
    S(t + 1, :) = S(t, :) + dS;
end

end