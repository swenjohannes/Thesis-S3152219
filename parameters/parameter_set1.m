%% Parameter set 1
S0 = 100;                           %Initial stock price
K = 80;                             %Strike price
L = 120;                            %Barrier level = L
v0 = 0.1;                           %Volatility
r = 0.05;                           %Risk-free interest rate
q = 0.02;                           %Divident yield
eta = 0.1;                          %Volatility of volatility
rho = 0.5;                          %Correlation coefficient
theta = 0.1;                        %Long-run mean of volatility   
kappa = 5;                          %Mean reversion speed
T = 1;                              %Time to maturity   
H = 0.1;                            %Roughness parameter
fprintf('Loaded parameter set 1')   %Succesful!