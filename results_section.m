%% Housekeeping
clear, clc;
close all
folder = fileparts(which(mfilename)); 
addpath(genpath(folder));              %load functions!

parameter_set2 %Load parameter set 1

%% Differences with respect to T and K

%Additional parameters
npath = 1e5; 
steps = 1200;
N = 100;
Nobs = 4;
K = 80:5:110;
T_vec = [8, 16, 32, 64, 92, 128, 250]/ 250;

%% MC simulation
[a2, b2, sim_res] = deal(NaN); %Start with empty values
[VT_c, VT_r] = deal(NaN([npath, length(T_vec)]));
for j = 1:length(T_vec)
    T = T_vec(j);
    [S_c, V_c] = heston_mc(S0, v0, rho, kappa, theta, T, r, q, eta, npath, steps);
    [S_r, V_r] = rough_heston_mc4(S0, v0, rho, kappa, theta, T, r, q, eta, H, npath, steps);
    
    %For later inspectation and to make cumulants
    VT_c(:, j) = sort(V_c(steps + 1, :));          
    VT_r(:, j) = sort(V_r(steps + 1, :));   
    a2(j) = VT_c(npath/1e3, j)           %Determine a2 as 0.001 tol!
    b2(j) = VT_r(npath - npath/1e3, j)   %Determine b2 as 0.001 tol!

    for k = 1:length(K)
        sim_res(j, k, 1) = barrier_prices_dm(S_c, K(k), L, Nobs, "uo", 1, r, T);
        sim_res(j, k, 2) = barrier_prices_dm(S_r, K(k), L, Nobs, "uo", 1, r, T);
    end
end
diff = sim_res(:, :, 1) - sim_res(:, :, 2);
[X, Y] = meshgrid(K, T_vec(1:6));
figure;
surf(X, Y, diff)

%Investigate a2 and b2
j = 3;
dx = 0.007;
x = (0:dx:1.5)';


T = T_vec(j);
C0 = eta ^ 2 * (1 - exp(- kappa * T)) / ( 4 * kappa); %Compute the constant
d = 4 * kappa * theta / eta ^ 2;                      %Degrees of freedom
lambda = 4 * kappa * exp(- kappa * T) * v0 ...        %Non-centrality ..
        / (eta ^ 2 * (1 - exp(- kappa * T)));         %parameter
ncx2 = ncx2pdf(x / C0,d,lambda);
%Create plot
figure;                                             
plot(x,ncx2,'b-','LineWidth',2); hold on;
[p,x] = hist(VT_c(:,j), x); plot(x,p/sum(p));             %Empirical pdf Heston
[p,x] = hist(VT_r(:,j), x); plot(x,p/sum(p));            %Empirical pdf rough Heston
legend('Non-central chi-squared distribution','Simulated Heston values','Simulated Rough Heston values')




%% Vary T and K: 2D-COS method
[S_r, V_r] = rough_heston_mc4(S0, v0, rho, kappa, theta, T, r, q, eta, H, npath, steps);
V_rh = sort(V_r(steps +1 , :));
a2 = V_rh(1e2);         
b2 = V_rh(1e5 - 1e2);  

N = 100;
cos_res = NaN;
for j = 1:length(T_vec)
    T = T_vec(j);
    for k = 1:length(K)
        cos_res(j, k, 1) = cos2d("h", S0, K(k), L, v0, r, q, eta, theta, rho, kappa, T, N, Nobs);
        cos_res(j, k, 2) = cos2d("rh", S0, K(k), L, v0, r, q, eta, theta, rho, kappa, T, N, Nobs, H,  a2, b2);
    end
end
cos_res
diff = cos_res(:, :, 1) - cos_res(:, :, 2);
[X, Y] = meshgrid(K, T_vec);
figure;
surf(X, Y, diff)
xlabel('Strike price') 
ylabel('Time to maturity') 
zlabel('Difference in option price') 
T_vec(3)

%% Vary the barrier level - keep K constant: 2D-COS method
parameter_set2;
N = 100;
cos_res = NaN;
L = 100:5:130;

for j = 1:length(T_vec)
    T = T_vec(j);
    for l = 1:length(L)
        cos_res(j, k, 1) = cos2d("h", S0, K, L(l), v0, r, q, eta, theta, rho, kappa, T, N, Nobs);
        cos_res(j, k, 2) = cos2d("rh", S0, K, L(l), v0, r, q, eta, theta, rho, kappa, T, N, Nobs, H,  0, 1.4);
    end
end
cos_res
diff = cos_res(:, :, 1) - cos_res(:, :, 2);
[X, Y] = meshgrid(L, T_vec);
figure;
surf(X, Y, diff)
xlabel('Barrier level') 
ylabel('Time to maturity') 
zlabel('Difference in option price') 

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

%MC simulation of Heston model
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

%H = 0.1;
%price = hrBarrier_mc(K,1000,1,1,1,T, 1, S0,r,q,v0,kappa,theta,eta,rho,H, npath, 250)
%[S_r, V_r, time_r] =  rough_heston_mc4(S0, v0, rho, kappa, theta, T, r, q, eta, H, npath, 250);
%price = barrier_prices_dm(S_r, K, 1000, 1, "uo", 1, r, T)
% cos2d2("rh", S0, K, 1000, v0, r, q, eta, theta, rho, kappa, 0.5, 120, 1, 0.1, 0, 1.4)

%cos2d2("rh", S0, K, 1000, v0, r, q, eta, theta, rho, kappa, 0.5, 120, 2, 0.1, 0, 1.4)

%Check vanilla price!
vanilla_prices(S_r, K, 1, r, T)
T =0.5;
H = 0.1;
cos2d2("rh", S0, K, 120, v0, r, q, eta, theta, rho, kappa, T, 100, 5, H, 0, 1.3, 1000)


%Determine a2 and b2 of rough Heston
V_rh = V_r(1201, : );                                
V_rh = sort(V_rh);
a2 = V_rh(1e2);         %0.0138 param set1
b2 = V_rh(1e5 - 1e2);   %0.2676 param set1




