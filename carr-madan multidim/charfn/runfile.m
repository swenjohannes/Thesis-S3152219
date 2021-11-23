%% Housekeeping
clear, clc;
close all
folder = fileparts(which(mfilename)); 
addpath(genpath(folder));              %load functions!

rd = 0.05;
rf = 0.00;
T = 1;
S0 = 100;
x0 = log(S0);
rho = -0.7;
sigma = 0.5; 
v0 = 0.04;
kappa = 1.5;
theta = 0.04;

n = 5; %number of observation points

H = 60;
K = 80; 

h = log(H);
k = log(K);

a = 1.2;
psi_N_ = @(v) psi_N(v, a, T, rho, sigma, kappa, theta, x0, v0, rd, rf);
psi_vec = @(v) vectorize_fun(v, psi_N_);

inner = integral(psi_vec, 0, 100);
exp(- a * k) / pi * real(inner)

phi_N(3,  T, rho, sigma, kappa, theta, x0, v0, rd, rf)
chfun_heston( rd, kappa, theta, v0, sigma, rho, T, 3, x0)


integral2(psi_N_, -100, 100, -100, 100)
