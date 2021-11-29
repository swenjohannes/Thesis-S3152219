function chi = chi_H(x1, x2, y, k, a, b)
pbma = pi / (b - a);
k_star = k * pbma;

%Arguments of cos/sin
c1 = k_star * (x1 - a);
c2 = k_star * (x2 - a);

%Compute Chi
if y == 0
    chi = x2 - x1;
else
    chi = 1 ./ (y ^ 2 + k_star .^ 2) ...
          .* (exp(x2 * y) .* (y * cos(c2) + k_star .* sin(c2)) ...
          -  exp(x1 * y) .* (y * cos(c1) + k_star .* sin(c1)));
end