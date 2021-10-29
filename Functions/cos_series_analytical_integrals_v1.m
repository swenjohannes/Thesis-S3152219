function [V_k_call, V_k_put] = cos_series_analytical_integrals_v1(a, b, k, K)

%{

 Based on code by Peter.Gruber@unisg.ch


 This function computes the Vk integrals for the Levy Processes and the Heston
 model. Vk coefficients enter the COS expansion for plain vanilla options. The following
 code is vectorised in k.
 Implementation follows Fang and Oosterlee (2008), Eq. (29) p. 7.


 Authors : Baldi Lanfranchi, Federico
         : La Cour, Peter

 Version : 1.0 (22.03.2019)


 [V_k_call, V_k_put] = cos_series_analytical_integrals_v1(a, b, k)


 Inputs : a             - Cosine argument (lower truncation bound) 
        : b             - Cosine argument (upper truncation bound)
        : k             - FFT Nx1 vector of k [0:N-1]


Outputs : V_k_call      - Nx1 vector of coefficients 
        : V_k_put       - Nx1 vector of coefficients 

%}


% Psi and Chi coeffiencts for the two option types
[psi_k_call, chi_k_call] = vanilla_cos_series_coef_v1(a, b, 0, b, k);
[psi_k_put, chi_k_put]   = vanilla_cos_series_coef_v1(a, b, a, 0, k);


% Vk coefficients 
V_k_call = 2 / (b - a) .* K' * (chi_k_call - psi_k_call);
V_k_put  = 2 / (b - a) .* K' * (- chi_k_put + psi_k_put);
