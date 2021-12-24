%% Parameter set 3
S0 = 100;                           %Initial stock price
K = 80;                             %Strike price
L = 120;                            %Barrier level = L
v0 = 0.3;                           %Volatility
r = 0.05;                           %Risk-free interest rate
q = 0.00;                           %Divident yield
eta = 0.5;                          %Volatility of volatility
rho = -0.7;                         %Correlation coefficient
theta = 0.3;                        %Long-run mean of volatility   
kappa = 3;                          %Mean reversion speed
T = 0.2;                            %Time to maturity   
H = 0.1;                            %Roughness parameter
fprintf('Loaded parameter set 3')   %Succesful!