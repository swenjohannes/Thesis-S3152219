%% Params
clear
clc

%% Parameters
N = 300;
k = 0:(N-1);

S0 = 100;
K = 80;
L = 110; %Barrier level = L

v0 = 0.2;
r = 0.00;
eta = 0.70;
rho = -0.7;
theta = 0.2;
kappa = 1.5;

x0 = log(S0 / K);
h = log(L / K);

T = 1;
Nobs = 2;
dt = T/ Nobs;

%Obtain a and b from heston cumulants
[c1, c2, ~]  = heston_cumulants_v1(r, kappa, theta, v0, eta, rho, T);
[a1, b1]= cos_truncation_range_v2(c1,c2,0,12);

%V_mc = sort(mc_V(v0, kappa, theta, T, eta, 1e6, 1024));
%a2 = V_mc(1000);
%b2 = V_mc(1e6 - 1000);

a2 = 0;
b2 = 1.5;

%% COS 2D
V(:, :, Nobs) = zeros(N);           %Create empty V dataframe
V(:, :, Nobs) = V_T(h, a1, b1, a2, b2, K, N); %Store V(T_M) at M


phi_Ap = phi_A(k, k, kappa, rho, eta, theta, r, dt, a1, b1, a2, b2);
phi_Am = phi_A(k, -k, kappa, rho, eta, theta, r, dt, a1, b1, a2, b2);

if Nobs > 1
    % %Obtain H positive and minus
    Hp = get_H(kappa, rho, eta, dt, 1, a1, b1, a2, b2, N);
    Hm = get_H(kappa, rho, eta, dt, -1, a1, b1, a2, b2, N);

    Mp = get_M(a1, h, a1, b1, N); %Obtain M

    for t = (Nobs - 1):-1:1
        Wp = 0.5 * exp(-r * dt) * phi_Ap .* V(:, :, t + 1); %W positive
        Wm = 0.5 * exp(-r * dt) * phi_Am .* V(:, :, t + 1); %W minus

        A = get_A(Hp, Hm, Wp, Wm, N);      %Obtain A
        A_ = A .* [0.5 repelem(1, N -1)]'; %Weight first row by a half!
        V(:, :, t) = real(Mp * A); %Compute the weighted sum in matrix form

    end
end
%Final step: compute V0

F = get_F(x0, v0, kappa, rho, eta, theta, r, dt, a1, b1, a2, b2, N);
V0 = exp(-r * dt) * weighted_sum_mat(0.5 * F .* V(:, :, 1))


