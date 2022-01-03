%% Housekeeping
clear, clc;
close all
folder = fileparts(which(mfilename)); 
addpath(genpath(folder));              %load functions!

parameter_set1 %Load parameter set 1
%parameter_set2 %Load parameter set 2

%Additional parameters
npath = 1e5; 
steps = 1200;
N = 100;

%MC simulation
[S_c, V_c] = heston_mc(S0, v0, rho, kappa, theta, T, r, q, eta, npath, steps);
[S_r, V_r] = rough_heston_mc4(S0, v0, rho, kappa, theta, T, r, q, eta, H, npath, steps);

%Simulation prices
barrier_prices_dm(S_c, K, L, Nobs, "uo", 1, r, T)
barrier_prices_dm(S_r, K, L, Nobs, "uo", 1, r, T)

%Cos method prices
cos2d("h", S0, K, L, v0, r, q, eta, theta, rho, kappa, T, N, Nobs)
%cos2d("rh", S0, K, L, v0, r, q, eta, theta, rho, kappa, T, N, Nobs, 0.6)
cos2d("rh", S0, K, L, v0, r, q, eta, theta, rho, kappa, T, N, Nobs, H, 0, 0.3)
cos2d2("rh", S0, K, L, v0, r, q, eta, theta, rho, kappa, T, N, Nobs, 1, 0, 0.5)
cos2d("h", S0, K, 1000, v0, r, q, eta, theta, rho, kappa, T, N, 1)
cos2d("rh", S0, K, 1000, v0, r, q, eta, theta, rho, kappa, T, N, 1, H)
vanilla_prices(S_c, K, 1, r, T)

%% Create results table 
parameter_set2;
K = 80;
T = 0.5;

%MC simulation of Heston model
npath = 1e5; 
steps = 1200;
Nobs = 2:5;
N = [40, 60, 80, 100];

[S_c, V_c, time_c] = heston_mc(S0, v0, rho, kappa, theta, T, r, q, eta, npath, steps);

results_c = NaN([(length(N)+ 1) length(Nobs)]);   %matrix to store results
results_c(1, :, 2) = time_c;          %store simulation time
for j = 1:length(Nobs)
    %Store simulation result in first column!
    results_c(1, j, 1) = barrier_prices_dm(S_c, K, L, Nobs(j), "uo", 1, r, T);
    for n = 1:length(N)
        [V0, time] = cos2d("h", S0, K, L, v0, r, q, eta, theta, rho, kappa, T, N(n), Nobs(j)); 
        results_c(n + 1, j, :) = [V0, time]; 
        fprintf('Completed inner step: %i', n)   
    end
    fprintf('Completed outer step: %i', j)   
end
results_c;

%MC simulation of Heston model
H = 0.1;
[S_r, V_r, time_r] = rough_heston_mc4(S0, v0, rho, kappa, theta, T, r, q, eta, H, npath, steps);

Nobs = 2:5;
N = [40, 60, 80, 100];
results_r = NaN([(length(N)+ 1) length(Nobs)]);   %matrix to store results

results_r(1, :, 2) = time_r;          %store simulation time
for j = 1:length(Nobs)
    %Store simulation result in first column!
    results_r(1, j, 1) = barrier_prices_dm(S_r, K, L, Nobs(j), "uo", 1, r, T); 
    for n = 1:length(N)
        [V0, time] = cos2d("rh", S0, K, L, v0, r, q, eta, theta, rho, kappa, T, N(n), Nobs(j), H, 0.001, 1.5, 300); 
        results_r(n + 1, j, :) = [V0, time]; 
        fprintf('Completed inner step: %i', n)   
    end
    fprintf('Completed outer step: %i', j)   
end
results_r;


results = horzcat(results_c, results_r)
values_results = results(:, :, 1); %Values
time_results = results(:, :, 2);    %times


cos2d("rh", S0, K, 1000, v0, r, q, eta, theta, rho, kappa, 0.3, N(n), 1, H, 0.001, 1.5, 500); 
vanilla_prices(S_r, K, 1, r, T)

%cos2d("rh", S0, K, L, v0, r, q, eta, theta, rho, kappa, T, N, j, H);

%% Rough heston test
%Truncation ranges
parameter_set2
T = 0.5;
[c1, c2, ~]  = heston_cumulants_v1(r, q, kappa, theta, v0, eta, rho, T);
[a1, b1] = cos_truncation_range_v2(c1,c2,0,12); %Obtain a and b from cumulants
[a2, b2] = a2_b2(eta, theta, kappa, T, v0, 1e-4);
N = 50;
Nobs = 2;
dt = T/Nobs;
k = 0:(N-1);        %Pre compute A and B2 parts of the characteristic equation!

a2 = 0, b2 = 0.3;

[phi_Ap, B2p]  = phi_A_B2_heston(k, k, kappa, rho, eta, theta, r, q, dt);
