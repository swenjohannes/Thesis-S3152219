%% Housekeeping
clear, clc;
close all
folder = fileparts(which(mfilename)); 
addpath(genpath(folder));              %load functions!

%% Parameters
%Model parameters
S0 = 100;           %initial price
V0 = 0.08;          %initial volatility
rho = -0.7;         
kappa = 1.5;        %mean reverting speed
theta = 0.06;       %mean of variance
T = 1;              %Time to maturity
r = 0.00;           %riskfree interest rate
nu = 0.7;           %vol param
H = 0.10;           %roughness

%Simulation parameters 
npath = 1e4;        %number of paths in the simulations
q = 11;             %to create multiple of 2 steps
N = 2^ q;           %number of simulation steps (multiple of 2 for FFT)

%% Simulation
%Classic heston - full Monte Carlo Simulation
name = 'Full MC simulation of classic Heston'
[S_c, V_c] = classic_heston_full_mc(S0, V0, rho, kappa, theta, T, r, nu, npath, N);
plot_SV(S_c, V_c, name);

%Euler scheme
name =  "Euler scheme simulation of rough heston"
[S_r, V_r] = rough_heston_euler_scheme(S0, V0, rho, kappa, theta, T, nu, H, npath, N);
plot_SV(S_r, V_r, name)

%% Compare Barier option prices
K = [80:1:120]; %strike prices
B = [80:1:120]; %barrier prices

%continiously monitored
close all
titles = ["Up and out", "Up and in", "Down and in", "Down and out"];
types = ["uo", "ui", "di", "do"];

fig_prices = figure; sgtitle("Prices of rough Heston versus Heston")
fig_abs_diff = figure; sgtitle("Absolute differences in prices of rough Heston versus Heston")
fig_rel_diff = figure; sgtitle("Relative differences in prices of rough Heston versus Heston")

for i = 1:length(types)
    prices_c = barrier_prices_cm(S_c, K, B, types(i), 1); %classic
    prices_r = barrier_prices_cm(S_r, K, B, types(i), 1); %rough
    
    
    abs_diff = prices_r - prices_c;
    rel_diff = abs_diff ./ prices_c;
    
    N = size(K, 2);
    M = size(B, 2);
    x = reshape(repelem(K, M), M, N)';
    y = reshape(repelem(B, N), N, M);

    figure(fig_abs_diff); %store these in 1 figure
    subplot(2, 2, i)
    surf(x, y, abs_diff,'edgecolor','none', 'FaceColor',[0.07 0.6 1], 'FaceAlpha',0.5); 
    zlabel('Absolute difference'); 
    ylabel('Barrier level');
    xlabel('Strike prices')
    title(titles(i))
    
    figure(fig_rel_diff); %store these in 1 figure
    subplot(2, 2, i)
    surf(x, y, rel_diff ,'edgecolor','none', 'FaceColor', [1 0 1] , 'FaceAlpha',0.5); 
    zlabel('Relative difference'); 
    ylabel('Barrier level');
    xlabel('Strike prices')
    title(titles(i))
    
    figure(fig_prices); %store these in 1 figure
    subplot(2, 2, i)
    surf(x, y, prices_r, 'edgecolor','none', 'FaceColor',[0.07 0.6 1], 'FaceAlpha',0.3);
    hold on;
    surf(x, y, prices_c, 'edgecolor','none', 'FaceColor',[1 0 1], 'FaceAlpha',0.3);
    zlabel('Price'); 
    ylabel('Barrier level');
    xlabel('Strike prices')
    title(titles(i));
end

%Checks