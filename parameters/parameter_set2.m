%% Parameter set 2
S0 = 100;                           %Initial stock price
K = 80;                             %Strike price
L = 120;                            %Barrier level = L
v0 = 0.2;                           %Volatility
r = 0.00;                           %Risk-free interest rate
q = 0.00;                           %Divident yield
eta = 0.5;                          %Volatility of volatility
rho = -0.7;                         %Correlation coefficient
theta = 0.22;                       %Long-run mean of volatility   
kappa = 1.7;                        %Mean reversion speed
T = 0.5;                             %Time to maturity   
H = 0.1;                            %Roughness parameter
disp('Loaded parameter set 2')   %Succesful!