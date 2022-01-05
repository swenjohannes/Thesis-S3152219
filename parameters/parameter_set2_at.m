if 0
    %% Parameter set 2
    S0 = 1.1225;                           %Initial stock price
    K = (1.08:0.01:1.16);
    x0 = log(S0./K);
    r = 0.00;                           %Risk-free interest rate
    q = 0.00;                           %Divident yield
    v0 = 0.006628295295875;                           %Volatility
    kappa = 1.513534749169168;                        %Mean reversion speed
    theta = 0.006682423753125;                       %Long-run mean of volatility   
    eta = 0.412024043065611;                          %Volatility of volatility
    rho = -0.216845113561275;                         %Correlation coefficient
    T = 0.5;                             %Time to maturity   
    H = 0.1;                            %Roughness parameter
    disp('Loaded parameter set 2 FX')   %Succesful!
else
    %% Parameter set 2
    S0 = 100;                           %Initial stock price
    K = 80;                             %Strike price
    L = 120;                            %Barrier level = L
    v0 = 0.2^2;                           %Volatility
    r = 0.00;                           %Risk-free interest rate
    q = 0.00;                           %Divident yield
    eta = 0.5;                          %Volatility of volatility
    rho = -0.7;                         %Correlation coefficient
    theta = 0.22;                       %Long-run mean of volatility   
    kappa = 1.7;                        %Mean reversion speed
    T = 0.5;                             %Time to maturity   
    H = 0.1;                            %Roughness parameter
    disp('Loaded parameter set 2')   %Succesful!
end