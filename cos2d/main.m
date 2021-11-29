%% Params
clear 
clc

N = 10;
k = 0:(N-1);

S0 = 100;
K = 80;
H = 1000;

v0 = 0.04;
r = 0.00;
eta = 0.7;
rho = -0.7;
theta = 0.04;
kappa = 1.5;

x0 = log(S0 / K);
h = log(H / K);

T = 1;
Nobs = 2;
dt = T/ Nobs;

[c1, c2, ~]  = heston_cumulants_v1(r, kappa, theta, v0, eta, rho, T);
[a, b]= cos_truncation_range_v2(c1,c2,0,12);

%Global constants used in all functions
%% COS 2D

%Create empty V dataframe
V(:, :, Nobs) = zeros(N);
V(:, :, Nobs) = V_T(h, a, b, K, N);

%Obtain H positive and minus
get_H_ = @(alpha) get_H(kappa, rho, eta, dt, alpha, a, b, N);
Hp = get_H_(1);
Hm = get_H_(-1);

%Obtain M
Mp = get_M(a, h, a, b, N);

for t = (Nobs - 1):-1:1
     %Obtain W positive and minus
    get_W_ = @(alpha) get_W(kappa, rho, eta, theta, r, dt, alpha, t, a, b, N, V);
    Wp = get_W_(1);
    Wm = get_W_(-1);
            
    %Obtain A
    A = get_A(Hp, Hm, Wp, Wm, N);

    for k1 = 1:N
        for k2 = 1:N          
            V(k1, k2, t) = real(weighted_sum(Mp(k1, :), A(:, k2)'));
        end
    end
end

%Final step: compute V0
F = F(x0, v0, kappa, rho, eta, theta, r, dt, a, b, N);

v0 = exp(-r * dt) * weighted_sum_mat(0.5 * F , V(:, :, 1))


