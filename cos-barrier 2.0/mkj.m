function m = mkj(k, j, x1, x2, a, b)
bma = b - a;
pbma = pi / bma;
integrand = @(x) exp( 1i * j * pbma * (x - a)) .* cos(k * pbma * (x - a));
m = 2 / bma * integral(integrand, x1, x2);
end