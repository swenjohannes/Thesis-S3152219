%Test swen
clear all
parameter_set2

T = 1;
[c1, c2, ~]  = heston_cumulants_v1(r, q, kappa, theta, v0, eta, rho, T);
[a1, b1] = cos_truncation_range_v2(c1,c2,0,12); %Obtain a and b from cumulants
[a2, b2] = a2_b2(eta, theta, kappa, T, v0, 1e-4);

N = 50;
Nobs = 1;
dt = T/Nobs;

w1 = (0:(N-1))' * pi / (b1 - a1);  %Turn w1 in w1* and make a column vector
w2 = (0:(N-1)) * pi / (b2 - a2);

H = 0.5;
[phi_Ap, B2p]  = phi_A_B2_heston(w1, w2, kappa, rho, eta, theta, r, q, dt);
[phi_Ap_r, B2p_r]  = phi_A_B2_rough_heston2(w1,w2, r -q,kappa,theta,eta,rho,H,dt, 100);


%% Create results table N/M
clear all
parameter_set2
T = 0.5;
%MC simulation of Heston model
npath = 1e5; 
steps = 1200;
Nobs = [2:5, 8, 12];
N = [60, 80, 100, 120];

[S_c, V_c, time_c] = heston_mc(S0, v0, rho, kappa, theta, T, r, q, eta, npath, steps);

results_c = NaN(length(Nobs), (length(N)+ 3));   %matrix to store results
results_c(:, 3, 2) = time_c;          %store simulation time
for j = 1:length(Nobs)
    %Store simulation result in first column!
    [price, lower, upper, ~] = barrier_prices_dm(S_c, K, L, Nobs(j), "uo", 1, r, T, 0.95); 
    results_c(j, 1, 1) = lower;
    results_c(j, 2, 1) = price;
    results_c(j, 3, 1) = upper;
    for n = 1:length(N)
        [V0, time] = cos2d2("h", S0, K, L, v0, r, q, eta, theta, rho, kappa, T, N(n), Nobs(j)); 
        results_c(j, n + 3, :) = [V0, time]; 
        fprintf('Completed inner step: %i', n)   
    end
    fprintf('Completed outer step: %i', j)   
end

%Nicer format table
rownames=arrayfun(@num2str, Nobs,'uni',0);
varnames_v = ['lower', 'mean', 'upper', arrayfun(@num2str, N,'uni',0)];
varnames_t = ['sim', arrayfun(@num2str, N,'uni',0)];
value_res_c =array2table(results_c(:,:,1), ...
            'VariableNames', varnames_v, 'RowNames', rownames)  
time_res_c =array2table(results_c(:,3:end,2), ...
            'VariableNames', varnames_t, 'RowNames', rownames)  

%MC simulation of rough Heston model
steps = 250;
H = 0.1;
T = 0.5;
[S_r, V_r, time_r] = rough_heston_mc4(S0, v0, rho, kappa, theta, T, r, q, eta, H, npath, steps);
results_r = NaN(length(Nobs), (length(N)+ 3));   %matrix to store results
results_r(:, 3, 2) = time_r;          %store simulation time
for j = 1:length(Nobs)
    %Store simulation result in first column!
     [price, lower, upper, ~] = barrier_prices_dm(S_r, K, L, Nobs(j), "uo", 1, r, T, 0.95); 
     results_r(j, 1, 1) = lower;
     results_r(j, 2, 1) = price;
     results_r(j, 3, 1) = upper;
    for n = 1:length(N)
        [V0, time] = cos2d2("rh", S0, K, L, v0, r, q, eta, theta, rho, kappa, T, N(n), Nobs(j), H, 0, 1.4, 500); 
        results_r(j, n + 3, :) = [V0, time]; 
        fprintf('Completed inner step: %i', n)   
    end
    fprintf('Completed outer step: %i', j)   
end

%Nicer format table
rownames=arrayfun(@num2str, Nobs,'uni',0);
varnames_v = ['lower', 'mean', 'upper', arrayfun(@num2str, N,'uni',0)];
varnames_t = ['sim', arrayfun(@num2str, N,'uni',0)];
value_res_r =array2table(results_r(:,:,1), ...
            'VariableNames', varnames_v, 'RowNames', rownames)  
time_res_r =array2table(results_r(:,3:end,2), ...
            'VariableNames', varnames_t, 'RowNames', rownames)  
%Bind results together
results = vertcat(results_c, results_r)
values_results = results(:, :, 1)  %Values
time_results = results(:, 3:end, 2)    %times








%% SOLVIE TEST
pbma1 = pi / (b1 - a1);
pbma2 = pi / (b2 - a2);
u1 = (0:50)'* pbma1;
u2 = (0:50)'* pbma2;

N = 50;
dt = T / N;
M = length(u1); 

% To define the Volterra integral equation:
c1 = - 0.5 *(u1 .^2 + 1i * u1);
beta = kappa - 1i * rho * eta * u1; 
c3 = 0.5*eta ^ 2;

[B2, phi_A] = deal(zeros([M, M]));
for j = u2
    j = 5;
    g = @(t) 1i * j * pbma2 * t.^-alpha / gamma(1-alpha);
    f = @(w) (c1 - beta .*w + c3*w.^2);
    [psi, Dalpha_psi, Dalpha_psi2] = SolveVIE2(f, g,  1i * j * pbma2, alpha,T, N, M); %solve VIE
%     close all;
    plot_imag(Dalpha_psi)
    plot_imag(Dalpha_psi2)
    % Integrate to get the characteristic function
%       B2(:, abs(j) + 1) = sum((Dalpha_psi(:, 1:end-1) + Dalpha_psi(:,2:end))/2, 2) *dt  ...
%                              + 1i * j * pbma2; %Don't know why this needs to be added.

    Dalpha_psi2(:, 1)
    B2(:, abs(j) + 1) = sum((Dalpha_psi2(:, 1:end-1) + Dalpha_psi2(:,2:end))/2, 2) *dt ...
                         + 1i * j * pbma2; %Don't know why this needs to be added.

    phi_A(:, abs(j) + 1) = exp(1i * u1 * (r - q) * T + ...
                    kappa * theta * sum((psi(:, 1:end-1) + psi(:, 2:end))/2, 2) *dt);
end
end