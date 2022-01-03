function [price, lower, upper, price_mean] = barrier_prices_dm(S, K, L, Nobs, type, call, r, T, cl, B)
%{
    Compute barrier prices from simulation results assuming 
    discrete time montioring.
    
    Input:      S       steps x paths simulation matrix 
                K       1 x 1 strike price
                L       1 x 1 barrier price
                Nobs     #of observation moments, will be spread out 
                        evenly over the time interval
                type    type of barrier option
                call    call = 1 or put = -1
                r       riskfree interest rate
                T       Time to maturity
                cl      (optional) confidence level    
                B       (optional) number of bootstrap draws  
    output:     prices  1 x 1  price of option
                lower   lower bootstrap confidence
                upper   upper bootstrap confidence
%}

npath = size(S, 2);
T_idx = size(S,1);           %final step index
t_obs = round(((1:Nobs) * (T_idx - 1) / Nobs) + 1); %index of observation moment
S = S(t_obs, :);    % filter observation moments

%Indicator function
switch type
    case 'uo'
        I = max(S, [], 1) <  L;
    case 'ui'
        I = max(S, [], 1) >= L;
    case 'do'
        I = min(S, [], 1) >= L;
    case 'di'
        I = min(S, [], 1) < L;
    otherwise
        error('Invalid type! Use 1 out ')
end

%Compute price
payoff = max(call * (S(Nobs, :) - K'), 0) .* I;  %regular call x indicator

%BOOTSTRAP CONFIDENCE IF NECCESSARY!
if (exist ( 'cl' , 'var' ))
    if( ~exist ( 'B' , 'var' )), B = 10000; end
    prices_B = NaN([1, B]); %empty vector to store bootstrapped prices
    for B_idx = 1:B
        sample_idx = randsample(npath, npath, true); %Draw a random sample of indexes
        prices_B(B_idx) = mean(payoff(:, sample_idx), 2) * exp(-r * T);
    end
    
    %Obtain bootstrap results
    prices_B = sort(prices_B);
    lower = prices_B(int16(B * (1- cl) + 1)); %lower 95% confidence interval!
    upper = prices_B(int16(B * cl + 1)); %upper 95% confidence interval!    
    price_mean  = mean(prices_B);
end  %Default value

price = mean(payoff, 2) * exp(-r * T);
end
