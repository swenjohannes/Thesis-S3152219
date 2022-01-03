function [x_opt, iv] = calibrate_heston(vol_dat, T, S0, K, r, q, N, x0)
%Calibrate: v0, eta, theta, rho, kappa
obj_fun = @(x) sse_loss(vol_dat, S0,  K, r, q, T, N, x(1),x(2), x(3), x(4), x(5));
lb = [0,0, 0, -1, 0];
ub = [1,3, 1, 0, 5];
A = [];
b = [];
Aeq = [];
beq = [];
x_opt = fmincon(obj_fun, x0,A,b,Aeq,beq,lb,ub);
[~, iv] = cos_1d("h", S0, x_opt(1), K, r, q, T, x_opt(2), x_opt(3), x_opt(4), x_opt(5), N); %Return IV's of optimal results
end

function SSE = sse_loss(vol_actual, S0, K, r, q, T, N, v0, eta, theta, rho, kappa)
    [~, iv] = cos_1d("h", S0, v0, K, r, q, T, eta, theta, rho, kappa, N);
    SSE = sum((vol_actual - iv).^2);
end