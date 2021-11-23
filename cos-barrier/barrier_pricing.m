%Model parameters
S0 = 100;           %initial price
V0 = 0.2;          %initial volatility
rho = -0.7;         
kappa = 1.5;        %mean reverting speed
theta = 0.04;       %mean of variance
T = 1;              %Time to maturity
r = 0.00;           %riskfree interest rate
nu = 0.7;           %vol param

N = 2^8;
k = [0:(N-1)];

H = 100;
K = 80; 
h = log(H/K);
x = log(S0/K);

M = 12; 
dt = T / M;

%% Initialization
[c1, c2, ~]  = bs_cumulants_v1(V0, r, T );
[a, b]= cos_truncation_range_v2(c1,c2,0,12);
bma = b- a;
pbma = pi / bma;

%compute Vk(tM)
V = zeros(N, M);
[psi_k, chi_k] = vanilla_cos_series_coef_v1(a, b, 0, h, k);
G = 2 / bma * K * (chi_k - psi_k);
V(:, M)= G;

%use up and out call
x1 = a; x2 = h; c = h; d = b;

%constuct ms and mc
j = (2 * N - 1):-1:0;
mc = mj(a, b, x1, x2, j);

j = [0:-1:(1 - N), nan, (N - 1):-1:1];
ms = mj(a, b, x1, x2, j);
ms(isnan(ms)) = 0; 

sgn = repmat([1, -1], 1, N);
d1 = fft(ms);
d2 = sgn .* fft(mc);

%char_fn = @(w) heston_char_fn(r, kappa, theta, V0, nu, rho, dt, w);
char_fn = @(w, t) heston_char_fn(r, kappa, theta, V0, nu, rho, t, w);

for m = M:-1:2 
    u = char_fn(k * pbma, T - (M - m) * dt) .* V(:, m)';
    u(1) = 0.5 * u(1);
    us = [u repelem(0, N)];     %pad N zero's to us
    
    Dus = fft(us);
    Msu = ifft(d1 .* Dus);
    Msu = Msu(1:N);
    Mcu = ifft(d2 .* Dus);
    Mcu = flip(Mcu(1:N));

    V(:, m - 1) = exp(-r * dt) / pi * imag(Msu + Mcu);
end

%calculate V0
inner = real(char_fn(k * pbma, dt) .* exp(1i * k * (x - a) * pbma));
inner(1) = 0.5 * inner(1);
price = exp(-r * dt) * sum(inner)


%Simulation price
[S_c, V_c] = classic_heston_full_mc_v2(S0, sqrt(V0), rho, kappa, theta, T, r, nu, npath, N);
barrier_prices_dm(S_c, 100, 110, 1, "uo", M)

%Check! 1 period
%prices = S_c(N + 1, :);
%idx = prices >= 110;
%prices(idx) = 0;
%mean(max(prices- 100, 0))
