function Ak = Ak_(tau, kappa, sigma, q_r, u_r, p_, n)
A_ = @(a, b) A(tau, a, b, kappa, sigma); %less input params

Ak = zeros(1, n + 1);  %initial value = 0
for k = 1:n             %compute A0 -> An = n +1 elems!
    a = q_r(k) + Ak(k);
    b = p_(sum(u_r(1:k)));
    Ak(k + 1) = A_(a, b);   %store A0 in A(1), A1 in A(2) etc
end
end

