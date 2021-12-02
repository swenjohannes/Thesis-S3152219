%% Params
clear 
clc

N = 100;
k = 0:(N-1);

S0 = 100;
K = 80;
L = 110; %Barrier level = L

v0 = 0.2;
r = 0.00;
eta = 0.7;
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
a2 = 0.001;
b2 = 3;

%% COS 2D
V(:, :, Nobs) = zeros(N);           %Create empty V dataframe
V(:, :, Nobs) = V_T(h, a1, b1, a2, b2, K, N); %Store V(T_M) at M
 
if Nobs > 1
    % %Obtain H positive and minus
    Hp = get_H(kappa, rho, eta, dt, 1, a1, b1, a2, b2, N);
    Hm = get_H(kappa, rho, eta, dt, -1, a1, b1, a2, b2, N);
    
    % %Obtain M
    Mp = get_M(a1, h, a1, b1, N);
    
    for t = (Nobs - 1):-1:1
         %Obtain W positive and minus
        Wp = get_W(kappa, rho, eta, theta, r, dt, 1, t,  a1, b1, a2, b2, N, V);
        Wm = get_W(kappa, rho, eta, theta, r, dt, -1, t,  a1, b1, a2, b2, N, V);
    
        %Obtain A
        A = get_A(Hp, Hm, Wp, Wm, N);
    
        for k1 = 1:N
            for k2 = 1:N          
                V(k1, k2, t) = real(weighted_sum(Mp(k1, :) .* A(:, k2)'));
            end
        end
    end
end
%Final step: compute V0

F = get_F(x0, v0, kappa, rho, eta, theta, r, dt, a1, b1, a2, b2, N);
V0 = exp(-r * dt) * weighted_sum_mat(0.5 * F .* V(:, :, 1))


