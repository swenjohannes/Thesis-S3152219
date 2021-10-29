%% Housekeeping
clear, clc;
close all
folder = fileparts(which(mfilename)); 
addpath(genpath(folder));              %load functions!

%Model parameters
S0 = 100;           %initial price
x = log(S0);        %log price
K = 85;             %strike price
V0 = 0.04;          %initial volatility
rho = -0.5;         %assume independence for the moment!
kappa = 1.5;        %mean reverting speed
theta = 0.06;       %mean of variance
T = 0.25;           %Time to maturity
r = 0.00;           %riskfree interest rate
nu = 0.7;           %vol param

int1 = @(w) real(exp(-i.*w*log(K)).*heston_char_fn(r, kappa, theta, V0, nu, rho, T, w-i, x)./ ...
    (i*w.*heston_char_fn(r, kappa, theta, V0, nu, rho, T, -i, x)));
int1 = integral(int1, 0,100); %numerical integration
pi1 = int1/pi+0.5;

 % Inner integral 2
int2 = @(w) real(exp(-i .* w * log(K)) .* heston_char_fn(r, kappa, theta, V0, nu, rho, T, w, x) ./ (i * w)) ;
int2 = integral(int2, 0,100); %numerical integration
pi2 = int2/pi+0.5; % final pi2

% 2nd step: calculate call value
heston_C = S0*pi1-exp(-r*T)*K*pi2

