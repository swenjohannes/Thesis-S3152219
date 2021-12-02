%% Params
clear 
clc

N = 20;
k = 0:(N-1);

S0 = 100;
K = 80;
L = 110; %Barrier level = L

v0 = 0.04;
r = 0.00;
eta = 0.7;
rho = -0.7;
theta = 0.04;
kappa = 1.5;

x0 = log(S0 / K);
h = log(L / K);

T = 1;
Nobs = 2;
dt = T/ Nobs;

%Obtain a and b from heston cumulants
[c1, c2, ~]  = heston_cumulants_v1(r, kappa, theta, v0, eta, rho, T);
[a, b]= cos_truncation_range_v2(c1,c2,0,12);

%Global constants used in all functions
%% COS 2D

V(:, :, Nobs) = zeros(N);           %Create empty V dataframe
V(:, :, Nobs) = V_T(h, a, b, 0, b, K, N); %Store V(T_M) at M
figure(1)
surf(V(:, :, Nobs))

a1 = a;
b1 = b;
a2 = 0;
b2 = b1;
ypoint = 1000;
kpoint = 20;
y1 = linspace(a1,b1,ypoint)';
dy1 = y1(2) - y1(1);
y2 = linspace(a2,b2,ypoint)';
dy2 = y2(2) - y2(1);
v_test = K*((exp(y1)-1).*(y1>0).*(y1<h))*ones(1,ypoint);
k1 = (0:kpoint-1)';
k2 = (0:kpoint-1)';
V_k1y2 = zeros(kpoint,ypoint);
for i = 1:ypoint
    V_k1y2(:,i) = cos(pi*k1*(y1-a1)'/(b1-a1))*v_test(:,i)*dy1;
end
V_k1k2 = zeros(kpoint,kpoint);
for j = 1:kpoint
    V_k1k2(j,:) = V_k1y2(j,:)*cos(pi*(y2-a2)*k2'/(b2-a2)) * 2/(b1-a1)*2/(b2-a2)*dy2;
end
figure(2)
surf(V_k1k2)
%figure(3)
%surf(V_k1y2)

%% Obtain H positive and minus
Hp = get_H(kappa, rho, eta, dt, 1, a, b, N);
Hm = get_H(kappa, rho, eta, dt, -1, a, b, N);

%Obtain M
Mp = get_M(a, h, a, b, N);

for t = (Nobs - 1):-1:1
     %Obtain W positive and minus
    Wp = get_W(kappa, rho, eta, theta, r, dt, 1, t, a, b, N, V);
    Wm = get_W(kappa, rho, eta, theta, r, dt, -1, t, a, b, N, V);

    %Obtain A
    A = get_A(Hp, Hm, Wp, Wm, N);

    for k1 = 1:N
        for k2 = 1:N          
            V(k1, k2, t) = real(weighted_sum(Mp(k1, :), A(:, k2)'));
        end
    end
end

%Final step: compute V0

F = get_F(x0, v0, kappa, rho, eta, theta, r, dt, a, b, N);
V0 = exp(-r * dt) * weighted_sum_mat(0.5 * F, V(:, :, 1))


