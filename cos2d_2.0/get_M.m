function M = get_M(u1, u2, a, b, N)
bma = a - b;
pbma = pi / bma;

M = zeros(N);
for k1 = 0:(N-1)
    for j1 = 0:(N-1)
    integrand = @(y) exp(1i * j1 * pbma * (y - a)) ...
                        .* cos(k1 * pbma * (y - a));
    M(k1 + 1, j1 + 1) = 2 / bma * integral(integrand, u1, u2);
    end
end

end
