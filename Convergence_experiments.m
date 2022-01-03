%% Test for #paths Heston
npath = 1e5; 
steps = 1200;

[S_c, ~, ~] = heston_mc(S0, v0, rho, kappa, theta, T, r, q, eta, npath, steps);

Nobs = 4;
t_obs = ((0:Nobs) *steps / Nobs) + 1; %index of observation moment
S_T = S_c(t_obs, :);    % filter observation moments

B = 1e4; %Bootstrap #
paths = [1e3, 2e3, 5e3, 8e3, 1e4, 3e4, 5e4, 8e4, 1e5]; %Points of graph
path_res = NaN([length(paths), 3]);

for path_idx = 1:length(paths)
    prices_B = NaN([1, B]); %empty vector to store bootstrapped prices
    for B_idx = 1:B
        sample_idx = randsample(npath, paths(path_idx), true); %Draw a random sample of indexes
        prices_B(B_idx) = barrier_prices_dm(S_T(:, sample_idx), K, L, Nobs, "uo", 1, r, T);
    end 
    %Obtain bootstrap results 
    prices_B = sort(prices_B);
    path_res(path_idx, 1) = prices_B((B / 100) * 5 + 1); %lower 95% confidence interval!
    path_res(path_idx, 2) = mean(prices_B);
    path_res(path_idx, 3) = prices_B((B / 100) * 95 + 1); %upper 95% confidence interval!
    fprintf('Completed path_idx %i \n', path_idx)   %Succesful!
end

figure;
plot(paths, path_res, '-o')
xlabel('Number of paths') 
ylabel('Option price') 
legend('Lower 95% confidence', 'Mean price', 'Upper 95% confidence')
path_res(6, 3) - path_res(6, 1)

%% Test for #paths rough Heston
[S_r, ~, ~] = rough_heston_mc4(S0, v0, rho, kappa, theta, T, r, q, eta, H, npath, steps);
t_obs = ((0:Nobs) *steps / Nobs) + 1; %index of observation moment
S_T = S_r(t_obs, :);    % filter observation moments

B = 10000; %Bootstrap #
paths = [1e3, 2e3, 5e3, 8e3, 1e4, 3e4, 5e4, 8e4, 1e5]; %Points of graph
path_res_r = NaN([length(paths), 3]);

for path_idx = 1:length(paths)
    prices_B = NaN([1, B]); %empty vector to store bootstrapped prices
    for B_idx = 1:B
        sample_idx = randsample(npath, paths(path_idx), true); %Draw a random sample of indexes
        prices_B(B_idx) = barrier_prices_dm(S_T(:, sample_idx), K, L, Nobs, "uo", 1, r, T);
    end 
    %Obtain bootstrap results 
    prices_B = sort(prices_B);
    path_res_r(path_idx, 1) = prices_B((B / 100) * 5 + 1); %lower 95% confidence interval!
    path_res_r(path_idx, 2) = mean(prices_B);
    path_res_r(path_idx, 3) = prices_B((B / 100) * 95 + 1); %upper 95% confidence interval!
    fprintf('Completed path_idx %i \n', path_idx)   %Succesful!
end
figure;
plot(paths, path_res)

%% Test for #steps
steps = [200, 400, 800, 1000, 1200, 1600, 2000];
npath = 5e4;
B = 10000; %Bootstrap #
Nobs = 4;

step_res_c = NaN([length(steps), 3]);
time_c = NaN;

for step_idx = 1:length(steps)
    [S_c, ~, time_c(step_idx)] = heston_mc(S0, v0, rho, kappa, theta, T, r, q, eta, npath, steps(step_idx));
     prices_B = NaN([1, B]); %empty vector to store bootstrapped prices
     t_obs = ((0:Nobs) *steps(step_idx) / Nobs) + 1; %index of observation moment
     S_T = S_c(t_obs, :);    % filter observation moments
     for B_idx = 1:B
         sample_idx = randsample(npath, npath, true); %Draw a random sample of indexes
         prices_B(B_idx) = barrier_prices_dm(S_T(:, sample_idx), K, L, Nobs, "uo", 1, r, T);
     end 
     %Obtain bootstrap results
     prices_B = sort(prices_B);
     step_res_c(step_idx, 1) = prices_B((B / 100) * 5 + 1); %lower 95% confidence interval!
     step_res_c(step_idx, 2) = mean(prices_B);
     step_res_c(step_idx, 3) = prices_B((B / 100) * 95 + 1); %upper 95% confidence interval!
     fprintf('Completed step_idx %i \n', step_idx)   %Succesful!
end 

step_res_r = NaN([length(steps), 3]);
time_r = NaN;
for step_idx = 1:length(steps)
    [S_r, ~, time_r(step_idx)] = rough_heston_mc4(S0, v0, rho, kappa, theta, T, r, q, eta, H, npath, steps(step_idx));
     prices_B = NaN([1, B]); %empty vector to store bootstrapped prices
     t_obs = ((0:Nobs) *steps(step_idx) / Nobs) + 1; %index of observation moment
     S_T = S_r(t_obs, :);    % filter observation moments
     for B_idx = 1:B
         sample_idx = randsample(npath, npath, true); %Draw a random sample of indexes
         prices_B(B_idx) = barrier_prices_dm(S_T(:, sample_idx), K, L, Nobs, "uo", 1, r, T);
     end 
     %Obtain bootstrap results
     prices_B = sort(prices_B);
     step_res_r(step_idx, 1) = prices_B((B / 100) * 5 + 1); %lower 95% confidence interval!
     step_res_r(step_idx, 2) = mean(prices_B);
     step_res_r(step_idx, 3) = prices_B((B / 100) * 95 + 1); %upper 95% confidence interval!
     fprintf('Completed step_idx %i \n', step_idx)   %Succesful!
end 


step_res_c
step_res_r

plot(step_res_c)
plot(step_res_r)



%% Second try steps of Heston!
steps = [200, 400, 800, 1000, 1200, 1600, 2000, 4000, 8000, 16000, 32000];
npath = 5e4;
B = 10000; %Bootstrap #
Nobs = 4;

step_res_c = NaN([length(steps), 3]);
time_c = NaN;
for step_idx = 1:length(steps)
     S_T = heston_mc_final_only(S0, v0, rho, kappa, theta, T, r, q, eta, npath, steps(step_idx), Nobs);
     
     prices_B = NaN([1, B]); %empty vector to store bootstrapped prices
     for B_idx = 1:B
         sample_idx = randsample(npath, npath, true); %Draw a random sample of indexes
         prices_B(B_idx) = barrier_prices_dm(S_T(:, sample_idx), K, L, Nobs, "uo", 1, r, T, false);
     end 
     %Obtain bootstrap results
     prices_B = sort(prices_B);
     step_res_c(step_idx, 1) = prices_B((B / 100) * 5 + 1); %lower 95% confidence interval!
     step_res_c(step_idx, 2) = mean(prices_B);
     step_res_c(step_idx, 3) = prices_B((B / 100) * 95 + 1); %upper 95% confidence interval!
     fprintf('Completed step_idx %i \n', step_idx)   %Succesful!
end 

plot(step_res_c)

















%% Fractional adams experiment
N = 50;
Nobs = 4;
int_steps = [80 120 160 200 250 500 1000];
res_frac = NaN([length(int_steps), 2]);
for j = 1:length(int_steps)
    [V0, time] = cos2d("rh", S0, K, 120, v0, r, q, eta, theta, rho, kappa, T, N, Nobs, H, int_steps(j));
    res_frac(j, : ) = [V0, time];
end
res_frac(:, 3) = abs(res_frac(:, 1) - res_frac(end, 1)); %Absolute error
res_frac = horzcat(int_steps', res_frac); 
res_frac



%% Error with respect to T
N = 100;
npath = 1e5; 
steps = 1200;
T_vec = [0.1, 0.11, 0.112, 0.115, 0.116, 0.117 0.12, 0.14, 0.2, 0.3];
abs_error = NaN([1, length(T_vec)]);
for T_idx = 1:length(T_vec)
    [S_c, ~, ~] = heston_mc(S0, v0, rho, kappa, theta, T_vec(T_idx), r, q, eta, npath, steps);
    v_sim = vanilla_prices(S_c, K, 1, r, T_vec(T_idx));
    v_cos = cos2d("h", S0, K, 1000, v0, r, q, eta, theta, rho, kappa, T_vec(T_idx), N, 1);
    abs_error(T_idx) = abs(v_cos - v_sim);
end
T_res = vertcat(T_vec, abs_error_1, abs_error);




%% Price differences with respect to T
T_vec = [8, 16, 32, 64, 92, 128]/ 250;

Nobs = 5;
L = 120;

K = 100;
[price_h, price_rh] = deal(NaN([1, length(T_vec)]));
for j = 1:length(T_vec)
    price_h(j) = cos2d("h", S0, K, L, v0, r, q, eta, theta, rho, kappa, T_vec(j), N, Nobs); 
    price_rh(j) = cos2d("rh", S0, K, L, v0, r, q, eta, theta, rho, kappa, T_vec(j), N, Nobs, H);
end
APE = abs(price_rh - price_h) ./ price_h * 100
plot(T_vec *250, APE)
plot(T_vec, price_rh, T_vec, price_h)

cos2d("rh", S0, 100, L, v0, r, q, eta, theta, rho, kappa, 1, N, Nobs, H, 400)
%MC simulation
npath = 1e5; 
steps = 1200;
j = 5;
[S_c, V_c] = heston_mc(S0, v0, rho, kappa, theta, T_vec(j), r, q, eta, npath, steps);
[S_r, V_r] = rough_heston_mc2(S0, v0, rho, kappa, theta, T_vec(j), r, q, eta, H, npath, steps);

barrier_prices_dm(S_c, K, L, Nobs, "uo", 1, r, T_vec(j))
barrier_prices_dm(S_r, K, L, Nobs, "uo", 1, r, T_vec(j))