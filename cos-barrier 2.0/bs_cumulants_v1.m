function [c1, c2, c4, w] = bs_cumulants_v1(u_0, mu, T )

%{
 
 Authors : Baldi Lanfranchi, Federico
         : La Cour, Peter

 Version : 1.0 (21.03.2019)

 [c1, c2, c4, w] = bs_cumulants_v1( u_0, mu, T )


 Inputs : mu            - log price drift rate
        : u_0           - mean (long run) volatility
        : T             - time to maturity


Outputs : c1             - first cumulant
        : c2             - second cumulant
        : c4             - fourth cumulant
        : w              - drift correction term

%}

c1 = mu * T;
c2 = u_0 * T;
c4 = 0;

% Drift correction term
w = 0;