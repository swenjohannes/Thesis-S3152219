function [prices] = barrier_prices_dm(S, K, B, obs, type, call, r, T)
%{
    Compute barrier prices from simulation results assuming 
    discrete time montioring.
    
    Input:      S       steps x paths simulation matrix 
                K       M x 1 vector of strike prices to calculate
                B       N x 1 vector of barrier prices to calculate
                obs     #of observation moments, will be spread out 
                        evenly over the time interval
                type    type of barrier option
                call    call = 1 or put = -1
                r       riskfree interest rate
                T       Time to maturity
    output:     prices  M x N table of prices 
%}

T_idx = size(S,1);           %final step index
M = size(K, 2);
N = size(B, 2);
prices = zeros(M, N);

t_obs = ((1:obs) * (T_idx - 1) / obs) + 1; %index of observation moment
S = S(t_obs, :);    % filter observation moments

%Indicator function
switch type
    case 'uo'
        I = max(S, [], 1) <  B';
    case 'ui'
        I = max(S, [], 1) >= B';
    case 'do'
        I = min(S, [], 1) >= B';
    case 'di'
        I = min(S, [], 1) < B';
    otherwise
        error('Invalid type! Use 1 out ')
end

%Compute prices
for i = 1:N
    payoff = max(call * (S(obs, :) - K'), 0) .* I(i, :);  %regular call x indicator
    prices(:, i) = mean(payoff, 2) * exp(-r * T);
end
