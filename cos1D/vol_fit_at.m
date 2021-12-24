% Read vol data
filename = 'vol_surface_2012.csv';
vol_tbl = readtable(filename);
K = vol_tbl{1,2:end}/10;
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

% Data
S0 = 1.1325; %Spot
r = 0;
q = 0;
v0 = 0.06^2;
kappa = 1.55;
theta = v0*2; 
eta = 0.17;
rho = -0.2;
x0 = [v0,kappa,theta,eta,rho];
pc = 1;
X = zeros(length(T),5);

for j = 1:length(T)
    Tc = T(j);
    vol = vol_dat(j,:);
    err_func = @(x) sum(sum((vol - get_iv(K,pc,Tc, S0,r,q,x)).^2));

    lb = [0.0001,0.0001, 0.0001, 0.0001, -1];
    ub = [1, 3, 10*v0, 1, 1];
    A = [];
    b = [];
    Aeq = [];
    beq = [];
    %get_iv(K,pc,Tc, S0,r,q,x0)
    x_opt = fmincon(err_func,x0,A,b,Aeq,beq,lb,ub);
    X(j,:) = x_opt;
end

%% 
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

%%
vol = vol_dat;
err_func = @(x) sum(sum((vol - get_iv(K,pc,T, S0,r,q,x)).^2));

lb = [0.0001,0.0001, 0.0001, 0.0001, -1];
ub = [1, 3, 10*v0, 1, 1];
A = [];
b = [];
Aeq = [];
beq = [];
%get_iv(K,pc,Tc, S0,r,q,x0)
x_opt = fmincon(err_func,x0,A,b,Aeq,beq,lb,ub);

%% 
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


