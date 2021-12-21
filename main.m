%% Housekeeping
clear, clc;
close all
folder = fileparts(which(mfilename)); 
addpath(genpath(folder));              %load functions!

parameter_set1 %Load parameter set

%Additional parameters
npath = 1e5; 
steps = 1200;
N = 100;

%MC simulation
[S_c, V_c] = heston_mc(S0, v0, rho, kappa, theta, T, r, q, eta, npath, steps);
[S_r, V_r] = rough_heston_mc2(S0, v0, rho, kappa, theta, T, r, q, eta, H, npath, steps);




barrier_prices_dm(S_c, K, L, Nobs, "uo", 1, r, T)
barrier_prices_dm(S_r, K, L, Nobs, "uo", 1, r, T)
%vanilla_prices(S_c, 80, 1, r, T)

Nobs = 2;
K = 80;
L = 120;
T = 1;

cos2d("h", S0, K, L, v0, r, q, eta, theta, rho, kappa, T, N, Nobs)
%cos2d("rh", S0, K, L, v0, r, q, eta, theta, rho, kappa, T, N, Nobs, 0.6)
cos2d("rh", S0, K, L, v0, r, q, eta, theta, rho, kappa, T, N, Nobs, 0.6)

%% Rough heston test
%Truncation ranges
[c1, c2, ~]  = heston_cumulants_v1(r, q, kappa, theta, v0, eta, rho, T);
[a1, b1] = cos_truncation_range_v2(c1,c2,0,12); %Obtain a and b from cumulants
[a2, b2] = a2_b2(eta, theta, kappa, T, v0, 1e-4);
dt = T/Nobs;
k = 0:(N-1);        %Pre compute A and B2 parts of the characteristic equation!

[phi_Ap, B2p]  = phi_A_B2_heston(k, k, kappa, rho, eta, theta, r, q, dt, a1, b1, a2, b2);
[phi_Am, B2m]  = phi_A_B2_heston(k, -k, kappa, rho, eta, theta, r, q, dt, a1, b1, a2, b2);

[phi_Ap_r, B2p_r]  = phi_A_B2_rough_heston(k, k, kappa, rho, eta, theta, r, q, a1, b1, a2, b2, 1, dt, 500);
[phi_Am_r, B2m_r]  = phi_A_B2_rough_heston(k, -k, kappa, rho, eta, theta, r, q, a1, b1, a2, b2, 1, dt, 500);

%Note: We need to add + j * pbma2 to the equation.. no idea why
