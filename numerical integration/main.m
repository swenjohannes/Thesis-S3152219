%params
S0 = 100;
K = 80;
sigma = 0.04;
T = 1;
r = 0.00;
a = 1.5;
k = log(K);
s0 = log(S0);

%% BSM model
%Functions to be passed 
ch_fun = @(w)chfun_norm(s0, sigma, r, T, w);
psi = @(v) Psi(v, ch_fun, r, T, a);

C_carr = cm_vannila(psi, k, a) %Numerical integration
[ctk, ku]  = cm_vannila_fft(psi);
ctk(870)
exp(ku(870))

%% Heston 
theta = 0.04;
kappa = 1.5;
nu = 0.7;
rho = -0.7;
ch_fun = @(w) chfun_heston( r, kappa, theta, sigma , nu, rho, T, w, s0);
psi = @(v) Psi(v, ch_fun, r, T, a);

C_hest = cm_vannila(psi, k, a) %Numerical integration

%% Simulation
%Classic heston - full Monte Carlo Simulation
[S_c, V_c] = classic_heston_full_mc_v2(S0, sigma, rho, kappa, theta, T, r, nu, 1e5, 1e3);
C_hest_sim = vanilla_prices(S_c, 80, 1)

