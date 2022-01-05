%% Housekeeping
folder = fileparts(which(mfilename)); 
addpath(genpath(folder));              %load functions!

%parameter_set2_at %Load parameter set
parameter_set2 %Load parameter set

%Additional parameters
npath = 2e4; 
steps = 200;
N = 100;
Nobs = 10;

%% MC simulation

% tic
% [S_c, V_c] = heston_mc(S0, v0, rho, kappa, theta, T, r, q, eta, npath, steps);
% prices = barrier_prices_dm(S_c, K, L, Nobs, "uo", 1, r, T);
% Fmc = mean(S_c(end,:));
% fprintf('Classic Heston MC: '),toc()

seed = rng;

tic
nobs = Nobs;
nstep = steps;
price = hsBarrier_mc(K,L,1,1,1,T,nobs, S0,r,q,v0,kappa,theta,eta,rho, npath,nstep,seed);
fprintf(['Classic Heston MC: ' num2str(price) '\t']), toc()


tic
nobs = Nobs;
nstep = steps;
price = hrBarrier_mc(K,L,1,1,1,T,nobs, S0,r,q,v0,kappa,theta,eta,rho,H, npath,nstep,seed);
F = S0*exp((r-q)*T);
fprintf(['Rough Heston MC: ' num2str(price) '\t']), toc()

tic
nobs = Nobs;
nstep = steps;
price = hrBarrier_mckr(K,L,1,1,1,T,nobs, S0,r,q,v0,kappa,theta,eta,rho,H, npath,nstep,seed);
F = S0*exp((r-q)*T);
fprintf(['Rough Heston MC reset kernel: ' num2str(price) '\t']), toc()


if 0

tic()
[S_r, V_r] = rough_heston_mc2(S0, v0, rho, kappa, theta, T, r, q, eta, H, npath, nstep);
prices = barrier_prices_dm(S_r, K, L, Nobs, "uo", 1, r, T)
mean(S_r(end,:))
toc()






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

[phi_Ap, B2p]  = phi_A_B2_heston(k, k, kappa, rho, eta, theta, r, q, dt);
[phi_Am, B2m]  = phi_A_B2_heston(k, -k, kappa, rho, eta, theta, r, q, dt);

[phi_Ap_r, B2p_r]  = phi_A_B2_rough_heston2(k, k, kappa, rho, eta, theta, r, q);
[phi_Am_r, B2m_r]  = phi_A_B2_rough_heston2(k, -k, kappa, rho, eta, theta, r, q);

%Note: We need to add + j * pbma2 to the equation.. no idea why
end