function [psi_k, chi_k] = vanilla_cos_series_coef_v1(a, b, c, d, k)

%{

 Based on code by Peter.Gruber@unisg.ch


 This function computes cosine series coefficients (Chi and Psi) for plain vanilla options. 
 Implementation follows Fang and Oosterlee (2008), Eqs. (22; 23) p. 6.


 Authors : Baldi Lanfranchi, Federico
         : La Cour, Peter

 Version : 1.0 (22.03.2019)


 phi_hes = heston_char_function_v1(mu, lambda, u_bar, u_0, eta, rho, omega, T)


 Inputs : a             - Cosine argument (lower truncation bound) 
        : b             - Cosine argument (upper truncation bound)
        : c             - Lower integration bound for Psi/Chi (Eqs. 20, 21)
        : d             - Upper integration bound for Psi/Chi (Eqs. 20, 21)
        : k             - FFT Nx1 vector of k [0:N-1]


Outputs : psi_k         - Nx1 vector of first cosine coefficients 
        : chi_k         - Nx1 vector of second cosine coefficients 

%}


% Preliminary variables
bma    = b - a;
uu     = k * pi / bma;

% Chi coefficient
chi_k = 1 ./ (1 + uu .^2 ) .* ...
      ( cos(uu * ( d - a))* exp(d) - cos(uu * (c - a))* exp(c) + uu .* sin(uu * (d - a)) * exp(d)- uu .* sin(uu*(c-a))* exp(c) );

% Psi coefficient
uu(1)  = 1;      % to avoid case differentiation (done 2 lines below)
psi_k    = 1 ./ uu .* ( sin(uu * ( d - a )) - sin( uu * ( c - a) ) );
psi_k(1) = d - c;
  
