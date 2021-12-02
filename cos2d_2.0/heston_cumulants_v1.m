function [c1, c2, w] = heston_cumulants_v1(mu, lambda, u_bar, u_0, eta, rho, T)

%{
 This code computes up to the 2nd cumulant of ln(St/K) for the Heston Model
 Equations are given in Fang and Oosterlee (2008), Table 11, p. 21

 Authors : Baldi Lanfranchi, Federico
         : La Cour, Peter

 Version : 1.0 (21.03.2019)

 [c1, c2, omega] = heston_cumulants(mu, lambda, u_bar, u_0, eta, rho, T)


 Inputs : mu            - log price drift rate
        : lambda        - speed of mean reversion
        : u_bar         - mean (long run) volatility
        : u_0           - initial volatility
        : eta           - volatility of the volatility (vol of vol)
        : rho           - correlation between Wiener processes (W1 and W2)
        : T             - time to maturity


Outputs : c1            - first cumulant (mean)
        : c2            - second cumulant (variance)
        : w             - drift correction term

%}


% First cumulant (mean)
c1 = mu * T + (1 - exp(- lambda * T)) * (u_bar - u_0) / (2 * lambda) - 0.5 * u_bar * T;

% Second cumulant (variance)
c2 = 0.125 / lambda ^ 3 * (eta * T * exp(- lambda * T) * (u_0 - u_bar) * (8 * lambda * rho - 4 * eta) + ...
     lambda * rho * eta * (1 - exp(- lambda * T)) * (16 * u_bar - 8 * u_0) + ...
     2 * u_bar * lambda * T * (- 4 * lambda * rho * eta + eta ^ 2 + 4 * lambda ^ 2) + ...
     eta ^ 2 * ((u_bar - 2 * u_0) * exp(-2 * lambda * T) + u_bar * (6 * exp(- lambda * T) - 7) + 2 * u_0) + ...
     8 * lambda ^ 2 * (u_0 - u_bar) * (1 - exp(-lambda * T)));
 
 
 % Drift correction term (0 for the Heston Model)
 w = 0;