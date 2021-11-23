function phi = phi_N(u, T, rho, sigma, kappa, theta, x0, v0, rd, rf)
n = length(u); %n equals the number of u's 

dt = T / n;
t = dt:dt:T;

h = x0 + (rd - rf) * t;
q = 1i * u * rho / sigma;
j = v0 + kappa * theta * t;

p_ = @(u) p(u, kappa, rho, sigma);
q_r = flip(q); %Flip to get reverse order to
u_r = flip(u); %make n-k+1 the first element

Ak = Ak_(dt, kappa, sigma, q_r, u_r, p_,  n);
Bk = Bk_(dt, kappa, sigma, theta, q_r, u_r, p_, n, Ak);

%Compute terms
term1 = dot(1i * u, h);
term2 = dot(q, j);
term3 = sum(Bk);
term4 = Ak(n + 1) * v0;

phi = exp(term1 - term2 + term3 + term4);
end