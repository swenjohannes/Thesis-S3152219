K = 80:5:100;
T = 0.2;
cos_heston(S0, v0, K, r, q, T, eta, theta, rho, kappa, 100)
cos_rough_heston(S0, v0, K, r, q, T, eta, theta, rho, kappa, 1000,  0.1)




%% Minimize SSE
filename = 'vol_surface_2012.csv';
vol_dat = csvread(filename, 1, 1) / 10;
S0 = 1.1325; %Spot
K = [10.949,10.992,11.034,11.077,11.119,11.162,11.204,11.246,11.289,11.331,11.374,11.416,11.459,11.501,11.544,11.586,11.628] / 10;
T = [1/52, 2/52, 1/12, 2/12, 3/12, 6/12, 9/12, 1, 2, 3, 5, 7, 10];

%vol_dat = vol_dat(5:13, :)
%T = T(5:13)

q = [0.05, 0.1, 0.2, 0.24, 0.52, 0.81, 1.14, 1.36, 1.43]/100
r = [-0.781, - 0.774, -0.767, -0.759, -0.723, -0.673, -0.549, -0.425, -0.274]/100
%Plot surface
v0 = 0.03;
theta = 0.1;

[price, iv] = deal(NaN(length(T), length(K)))
for t = 1:length(T)
    [price(t, :), iv(t, :)] = cos_price(S0, v0, K, r, q, T(t), eta, theta, rho, kappa, 200)
end


%% Fit surface

%fit 5 params
%q = 0.0001;
%r = 0;
v0 = 0.2;
rho = -0.7;
obj_fun = @(x) sse(vol_dat, S0, x(1), K, r, q, T, x(2), x(3), x(4), x(5));
x0 = [v0, eta, theta, rho, kappa];
obj_fun(x0)

lb = [0,0, 0, -1, 0];
ub = [1,1, 1, 0, 5];
A = [];
b = [];
Aeq = [];
beq = [];
x_opt = fmincon(obj_fun, x0,A,b,Aeq,beq,lb,ub)
obj_fun(x_opt)
iv = get_iv_surf(S0, x_opt(1), K, r, q, T, x_opt(2), x_opt(3), x_opt(4), x_opt(5))


v0 = 0.2;
theta = 0.4;
r = 0.05;
q = 0.02;
kappa = 2;
rho  = -0.7;
eta = 0.4;
iv = get_iv_surf(S0, v0, K, r, q, T, eta, theta, rho, kappa)
close all
figure;
[X, Y] = meshgrid(K, T);
surf(X, Y, vol_dat); hold on
surf(X, Y, iv)

surf(vol_dat); hold on
surf(iv)
%figure;
%[X, Y] = meshgrid(K, T)
%surf(X, Y, iv)

 blsimpv(4722.35,4075,-0.02,32/250,650.15)


 K = [0.9:0.02:1.3]
 T = [1/52, 2/52, 1/12, 2/12, 3/12, 6/12, 9/12, 1, 2]
%% Vol surface compare
v0 = 0.4;
theta = 0.4;
r = 0.05;
q = 0.0;
kappa = 1;
rho  = -0.4;
eta = 0.5;
N =200;
H = 0.1;

iv_h = get_iv_surf("h", S0, v0, K, r, q, T, eta, theta, rho, kappa, N);
iv_rh = get_iv_surf("rh", S0, v0, K, r, q, T, eta, theta, rho, kappa, N, H);

close all
figure;
[X, Y] = meshgrid(K, T);
surf(X, Y, iv_h); hold on
surf(X, Y, iv_rh);
%surf(X, Y, vol_dat)
