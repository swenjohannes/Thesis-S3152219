function Bk = Bk_(tau, kappa, sigma, theta, q_r, u_r, p_, n, Ak)
B_ = @(a, b) B(tau, a, b, kappa, sigma, theta); %less input params

Bk = zeros(1, n); %initial value = 0
for k = 1:n
    a = q_r(k) + Ak(k); % note: Ak-1 is stored at location k!
    b = p_(sum(u_r(1:k)));
    Bk(k) = B_(a, b); %store B1 in B(1), B2 in B(2) etc
end

end