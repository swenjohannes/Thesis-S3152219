%% Housekeeping
clear, clc;
close all
folder = fileparts(which(mfilename)); 
addpath(genpath(folder));              %load functions!

parameter_set1 %Load parameter set

%Additional parameters
npath = 4e5; 
steps = 1200;
N = 300;

%MC simulation
[S_c, V_c] = heston_mc(S0, v0, rho, kappa, theta, T, r, q, eta, npath, steps);
[S_r, V_r] = rough_heston_mc(S0, v0, rho, kappa, theta, T, eta, H, npath, steps);

N_obs = 5;
cos2d_hest(S0, K, L, v0, r, q, eta, theta, rho, kappa, T, N, N_obs)
barrier_prices_dm(S_c, 80, 120, N_obs, "uo", 1, r, T)
%vanilla_prices(S_c, 80, 1, r, T) 

