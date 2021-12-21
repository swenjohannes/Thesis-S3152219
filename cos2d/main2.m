%% Params
clear 
clc
close all

%% Parameter set 1
N = 100;

S0 = 100;
K = 80;
L = 120; %Barrier level = L

v0 = 0.2;
r = 0.05;
q = 0.02;
eta = 0.7;
rho = 0.5;
theta = 0.1;
kappa = 1.5;

x0 = log(S0 / K);
h = log(L / K);

T = 1;
Nobs = 5;
dt = T/ Nobs;

%Obtain a and b from heston cumulants
[c1, c2, ~]  = heston_cumulants_v1(r, q, kappa, theta, v0, eta, rho, T);
[a1, b1] = cos_truncation_range_v2(c1,c2,0,12);

% V_mc = sort(mc_V(v0, kappa, theta, T, eta, 1e6, 1024));
% a2 = V_mc(1e3);
% b2 = V_mc(1e6 - 1e3);

a2 = 0;
b2 = 0.5;

%% COS 2D
V(:, :, Nobs) = zeros(N);           %Create empty V dataframe
V(:, :, Nobs) = V_T(h, a1, b1, a2, b2, K, N); %Store V(T_M) at M

% VT_an = V(:, :, Nobs);
% figure(1)
% surf(VT_an)


ypoint = 100;
% kpoint = N;
y1 = linspace(a1,b1,ypoint)';
y2 = linspace(a2,b2,ypoint)';
% V_y1y2 = K*((exp(y1)-1).*(y1>0).*(y1<h))*ones(1,ypoint);
% k1 = (0:kpoint-1)';
% k2 = (0:kpoint-1)';
% 
% V_k1k2 = cft2d(V_y1y2,kpoint,kpoint);
% figure(2)
% surf(V_k1k2)
% 
% v_y1y2 = icft2d(V_k1k2,ypoint,ypoint);
% 
% figure(3)
% surf(V_y1y2)
% figure(4)
% surf(v_y1y2)
% 

v_y1y2 = icft2d(V(:, :, Nobs),ypoint,ypoint);
figure();
surf(v_y1y2);

%Pre compute A and B2 parts of the characteristic equation!
k = 0:(N-1);
[phi_Ap, B2p]  = phi_A_B2(k, k, kappa, rho, eta, theta, r, q, dt, a1, b1, a2, b2);
[phi_Am, B2m]  = phi_A_B2(k, -k, kappa, rho, eta, theta, r, q, dt, a1, b1, a2, b2);

if Nobs > 1
    % %Obtain H positive and minus
    Hp = get_H(1, B2p, a2, b2, N);
    Hm = get_H(-1, B2p, a2, b2, N);

    Mp = get_M(a1, h, a1, b1, N); %Obtain M

    for t = (Nobs - 1):-1:1
        Wp = 0.5 * exp(-r * dt) * phi_Ap .* V(:, :, t + 1); %W positive
        Wm = 0.5 * exp(-r * dt) * phi_Am .* V(:, :, t + 1); %W minus

        A = get_A(Hp, Hm, Wp, Wm, N);      %Obtain A
        A(1, :) = A(1, :) * 0.5; %Weight first row by a half!
        V(:, :, t) = real(Mp * A); %Compute the weighted sum in matrix form
    end
end

%Plot of last Vk1k2's 
h= figure();
subplot(1, 2, 1); surf(V(:, :, t));
v_y1y2 = icft2d(V(:, :, t),ypoint,ypoint);
subplot(1, 2, 2); surf(v_y1y2);

%Final step: compute V0
F = get_F2(x0, v0, phi_Ap, phi_Am, B2p, B2m, a1, b1, a2, b2, N);
V0 = exp(-r * dt) * weighted_sum_mat(F .* V(:, :, 1))
 
%Calculate V0 as the weigthed sum (2D)
F(1, :) = F(1, :) * 0.5; %weight first row a by a half
F(:, 1) = F(:, 1) * 0.5; %weigth first col by a half
V0 = exp(-r * dt)  * sum(F .* V(:, : , 1), 'all')




