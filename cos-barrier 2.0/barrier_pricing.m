%% Params
N = 128;
k = 0:(N-1);

S0 = 100;
K = 80;
H = 110;

sigma = 0.2;
r = 0.00;
eta = 0.7;
rho = -0.7;
theta = 0.2;
kappa = 1.5;

x = log(S0 / K);
h = log(H / K);

T = 1;
M = 16;
dt = T/ M;

%% Compute a/b 
[c1, c2, ~]  = heston_cumulants_v1(r, kappa, theta, sigma, eta, rho, T);
[a, b]= cos_truncation_range_v2(c1,c2,0,12);
bma = b - a;
pbma = pi /bma;

x1 = a; x2 = h; c = h; d = b;
Mkj = mkj_mat(N, x1, x2, a, b);

%% Compute Vk(tM) according to formula 32
V = zeros(N, M);
[psi_k, chi_k] = vanilla_cos_series_coef_v1(a, b, 0, h, k);

V(:, M) = 2 / bma * K * (chi_k - psi_k);

%char_fn = @(j, t, sigma) heston_char_fn(r, kappa, theta, sigma, eta, rho, t, j *pbma);
char_fn = @(j, t) chfun_norm(sigma, r, t, j * pbma);

vt = ones(1, M) * sigma;
for m = 2:M  
    vt(m) = vt(m - 1) + kappa * (theta - vt(m -1)) * dt;
end
vt; 

%Recursively compute Vk(tm) from M - 1 to 1
for m = (M - 1):-1:1
    %phi = char_fn(k, dt, vt(m + 1))';
    phi = char_fn(k, dt)';
    for n = 1:N
        f = phi .* V(:, m + 1) .* Mkj(n, :)'; %multiply elementwise
        f(1) = 0.5 * f(1); %weight first element half
        V(n, m ) = exp(-r * dt) * real (sum(f));
    end
end

%Compute V0:
%phi = char_fn(k, dt, vt(1))';
phi = char_fn(k, dt)';
f = phi .* exp(1i * k * pbma * (x - a))' .* V(:, 1);
f(1) = 0.5 * f(1);
v0 = exp(- r * dt) * sum(real(f))

%Regular call:
% [psi_k, chi_k] = vanilla_cos_series_coef_v1(a, b, 0, b, k);
% UK = 2 / bma  * (chi_k - psi_k);
% phi = char_fn(k, T)';
% q = log(S0 / K);
% f = phi .* exp(1i * k * pbma * (q - a))';
% f(1) = 0.5 * f(1);
% v0_c = K * exp(- r * dt) * sum(real(f'  .* UK)) 



%% Simulation

%Simulation parameters 
npath = 1e5;        %number of paths in the simulations
q = 11;             %to create multiple of 2 steps
nsteps = 2^ q;           %number of simulation steps (multiple of 2 for FFT)

%BS
[S_bs] = bs_mc(S0, sigma, r, T, npath, nsteps);
barrier_prices_dm(S_bs, K, H, M, "uo", 1)
vanilla_prices(S_bs, 80, 1) %To check H -> inf and M = 1!
 
%Heston
[S_c, V_c] = classic_heston_full_mc_v2(S0, sigma, rho, kappa, theta, T, r, eta, npath, nsteps);
barrier_prices_dm(S_c, K, H, M, "uo", 1)
vanilla_prices(S_c, 80, 1) %To check H -> inf and M = 1!

