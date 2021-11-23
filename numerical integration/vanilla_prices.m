function [prices] = vanilla_prices(S, K, call)
%{
    Compute barrier prices from simulation results assuming 
    continious montioring.
    
    Input:      S       steps x paths simulation matrix 
                K       M x 1 vector of strike prices to calculate
                call    call = 1 or put = -1
    output:     prices  M vector of prices 
%}

T = size(S,1);           %final step index
M = size(K, 2);
prices = zeros(M, 1);

%Compute prices
payoff = max(call * (S(T, :) - K'), 0);  %regular call x indicator
prices = mean(payoff, 2);
