% Read vol data
filename = 'aex_vol_2412.csv';
vol_tbl = readtable(filename);
K = vol_tbl{1,2:end};
Tvec = str2tenor(vol_tbl{2:end,1});
T = Tvec *[1,1/12,1/365,0,0,0]';
vol_dat = vol_tbl{2:end,2:end}/100;

figure(1)
[X, Y] = meshgrid(K, T);
%scatter3(X(:), Y(:), vol_dat(:));
surf(X, Y, vol_dat);
xlabel('Strike')
ylabel('Maturity')
zlabel('Volatility')

%% Filter maturities between 2W and 3W and strikes between 0.75 and 1.25
T = T(1:13);
K = K(4:14);
vol_dat = vol_dat(1:13, 4:14)

figure(1)
[X, Y] = meshgrid(K, T);
%scatter3(X(:), Y(:), vol_dat(:));
surf(X, Y, vol_dat);
xlabel('Strike')
ylabel('Maturity')
zlabel('Volatility')

% Data
S0 = 1; %Spot
r = 0;
q = 0;
v0 = 0.2^2;
kappa = 5;
theta = v0*2; 
eta = 1.5;
rho = -0.6;
x0 = [v0,kappa,theta,eta,rho];
pc = 1;
X = zeros(length(T),5);

for j = 1:length(T)
    Tc = T(j);
    vol = vol_dat(j,:);
    err_func = @(x) sum(sum((vol - get_iv(K,pc,Tc, S0,r,q,x)).^2));
    lb = [0.0001,0.0001, 0.0001, 0.0001, -1];
    ub = [1, 3, 10*v0, 3, 1];
    A = [];
    b = [];
    Aeq = [];
    beq = [];
    %get_iv(K,pc,Tc, S0,r,q,x0)
    x_opt = fmincon(err_func,x0,A,b,Aeq,beq,lb,ub);
    X(j,:) = x_opt;
end
Xc = X;

%%  Classic Heston per maturity
figure(2)
%hold on
col = hsv(length(T));
for j = 1:length(T)
    plot(K,vol_dat(j,:),'o', K,get_iv(K,pc,T(j), S0,r,q,X(j,:)),'Color',col(j,:))
    hold all
end
xlabel('Strike')
ylabel('Vol')
title('Classic Heston per maturity')
hold off

%% Global fit
x0 = Xc(6, :) % x0 = [0.0590,    4.9882,    0.0394,    1.3953,   -0.7098]
vol = vol_dat;
err_func = @(x) sum(sum((vol - get_iv(K,pc,T, S0,r,q,x)).^2));

lb = [0.0001,0.0001, 0.0001, 0.0001, -1];
ub = [1, 5, 10*v0, 3, 1];
A = [];
b = [];
Aeq = [];
beq = [];
%get_iv(K,pc,Tc, S0,r,q,x0)
x_opt = fmincon(err_func,x0,A,b,Aeq,beq,lb,ub);

%% Plot global fit
figure(3)
%hold on
col = hsv(length(T));
for j = 1:length(T)
    plot(K,vol_dat(j,:),'o', K,get_iv(K,pc,T(j), S0,r,q,x_opt),'Color',col(j,:))
    hold all
end
xlabel('Strike')
ylabel('Vol')
title('Classic Heston global')
hold off

x_opt = [0.0323    2.7587    0.0252    0.7025   -0.6400]
figure;
[X, Y] = meshgrid(K, T);
surf(X, Y, vol_dat); hold on;
surf(X, Y, get_iv(K,pc,T, S0,r,q,x_opt));

%% CONTINUE HERE! FIT ROUGH HESTON!

for j = 1:length(T)
    get_iv_hr(K,pc,T(j), S0,r,q, x)
    disp(j)
end

%% Rough Heston
% Data
x0 = [0.0323    2.7587    0.0252    0.7025   -0.6400]

err_func = @(x) sum(sum((vol - get_iv_hr(K,pc,T, S0,r,q,x)).^2));

lb = [0.0001,0.0001, 0.0001, 0.0001, -1, 0.1];
ub = [1, 5, 10*v0, 3, 1, 0.2];
A = [];
b = [];
Aeq = [];
beq = [];

[x, fval] = fmincon(err_func,x0,A,b,Aeq,beq,lb,ub,@confuneq);


x= [ 0.0416    3    0.0462    0.6   -0.7    0.1]
%get_iv_hr(K,pc,T, S0,r,q,x)

figure;
[X, Y] = meshgrid(K, T);
%surf(X, Y, vol_dat); hold on;
surf(X, Y, get_iv_hr(K,pc,T, S0,r,q,x))





%%  Rough Heston per maturity
figure(2) 
col = hsv(length(T));
for j = 1:length(T)
    j =1
    plot(K,vol_dat(j,:),'o', K,get_iv_hr(K,pc,T(j), S0,r,q,X(j,:)),'Color',col(j,:))
    hold all
end
xlabel('Strike')
ylabel('Vol')
title('Rough Heston per maturity')
hold off


%% fit a single maturity

% Data
S0 = 1; %Spot
r = 0;
q = 0;
v0 = 0.18^2;
kappa = 1;
theta = v0* 3; 
eta = 1.5;
rho = -0.5;
pc = 1;

x0 = [v0,kappa,theta,eta,rho];

j = 15;
Tc = T(j);
iv = get_iv(K,pc,Tc, S0,r,q,x0);

figure;
plot(K, get_iv_hr(K,pc,Tc, S0,r,q,[x0, 0.2])); hold on;
plot(K, vol_dat(j, :)); %plot smile at T = j
title("Smile at time j")
xlabel("Moneyness")
ylabel("Implied vol")
legend('iv', 'aex iv')

err_func = @(x) sum(sum((vol - get_iv(K,pc,Tc, S0,r,q,x)).^2));
x0 = [v0,kappa,theta,eta,rho];
lb = [0.0001,0.0001, 0.0001, 0.0001, -1];
ub = [1, 3, 10*v0, 4, 1];
A = [];
b = [];
Aeq = [];
beq = [];
x_opt = fmincon(err_func,x0,A,b,Aeq,beq,lb,ub);





