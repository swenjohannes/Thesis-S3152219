function [prices] = barrier_prices_cm(S, K, B, type, call)
%{
    Compute barrier prices from simulation results using
    continious montioring.
    
    Input:      S       steps x paths simulation matrix 
                K       M x 1 vector of strike prices to calculate
                B       N x 1 vector of barrier prices to calculate
                type    type of barrier option
                call    call = 1 or put = -1
    output:     prices  M x N table of prices 
%}

T = size(S,1);           %final step index
M = size(K, 2);
N = size(B, 2);
prices = zeros(M, N);

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
    payoff = max(call * (S(T, :) - K'), 0) .* I(i, :);  %regular call x indicator
    prices(:, i) = mean(payoff, 2);
end
