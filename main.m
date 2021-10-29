%% Housekeeping
clear, clc;
close all
folder = fileparts(which(mfilename)); 
addpath(genpath(folder));              %load functions!

%% Parameters
%Model parameters
S0 = 100;           %initial price
V0 = 0.04;          %initial volatility
K = 85;             %strike price
rho = -0.5;         %assume independence for the moment!
kappa = 1.5;          %mean reverting speed
theta = 0.06;       %mean of variance
T = 0.25;              %Time to maturity
r = 0.00;           %riskfree interest rate
nu = 0.7;           %vol param
H = 0.15;            %roughness

%Simulation parameters 
npath = 1e4;        %number of paths in the simulations
q = 10;             %to create multiple of 2 steps
N = 2^ q;           %number of simulation steps (multiple of 2 for FFT)

%% Simulation
%Classic heston - full Monte Carlo Simulation
name = 'Full MC simulation of classic Heston'
[S, V] = classic_heston_full_mc(S0, V0, rho, kappa, theta, T, r, nu, npath, N);
plot_SV(S, V, name);

%Classic heston - full Monte Carlo Simulation
name = 'Full MC simulation of rough Heston'
[S, V] = rough_heston_full_mc(S0, V0, rho, kappa, theta, T, r, nu, H, npath, N);
plot_SV(S, V, name)


%This one tends to explode.. don't know why 
name = "Approx MC simulation of rough heston"
[S, V] = rough_heston_approx_mc(S0, V0, rho, kappa, theta, T, r, nu, H, npath, N);
plot_SV(S, V, name)

%Euler scheme
name =  "Euler scheme simulation of rough heston"
[S, V] = rough_heston_euler_scheme(S0, V0, rho, kappa, theta, T, nu, H, npath, N);
plot_SV(S, V, name)

%Call value!
call_value = mean(max(S(N+1, :) - 85, 0) ) 

%% Semi analytic solution
