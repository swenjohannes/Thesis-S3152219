function [c1, c2, w] = heston_cumulants_v1(r, q, kappa, theta, v0, eta, rho, tau)

%{
 Description: This code computes up to the 2nd cumulant of ln(St/K) for the Heston Model
               Equations are given in Fang and Oosterlee (2008), Table 11, p. 21
               Original code from : Baldi Lanfranchi, Federico, La Cour, Peter


 [c1, c2, omega] = heston_cumulants(r, kappa, theta, v0, eta, rho, tau)


 Inputs : r            - riskfree interest rate
        : q            - dividend yield
        : kappa        - speed of mean reversion
        : theta        - mean (long run) volatility
        : v0           - initial volatility
        : eta          - volatility of the volatility (vol of vol)
        : rho          - correlation between Wiener processes (W1 and W2)
        : tau          - time to maturity


Outputs : c1            - first cumulant (mean)
        : c2            - second cumulant (variance)
        : w             - drift correction term
%}


% First cumulant (mean)
c1 = (r - q) * tau + (1 - exp(- kappa * tau)) * (theta - v0) / (2 * kappa) - 0.5 * theta * tau;

% Second cumulant (variance)
c2 = 0.125 / kappa ^ 3 * (eta * tau * exp(- kappa * tau) * (v0 - theta) * (8 * kappa * rho - 4 * eta) + ...
     kappa * rho * eta * (1 - exp(- kappa * tau)) * (16 * theta - 8 * v0) + ...
     2 * theta * kappa * tau * (- 4 * kappa * rho * eta + eta ^ 2 + 4 * kappa ^ 2) + ...
     eta ^ 2 * ((theta - 2 * v0) * exp(-2 * kappa * tau) + theta * (6 * exp(- kappa * tau) - 7) + 2 * v0) + ...
     8 * kappa ^ 2 * (v0 - theta) * (1 - exp(-kappa * tau)));
 
 
 % Drift correction term (0 for the Heston Model)
 w = 0;