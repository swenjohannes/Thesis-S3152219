function chi = chi_H(x1, x2, B, a, b, N)

pbma = pi / (b - a);
k = zeros(1); %Create 1x1 element
k(1, 1, N) = 0; %Add 3D dimension
k(1, 1, :) =  (0:(N-1)) * pbma; %Give proper values

%Arguments of cos/sin
c1 = k * (x1 - a);
c2 = k * (x2 - a);

chi = 1 ./ (B .^2 + k .^ 2) ...
    .* ( exp(x2 * B) .* (B .* cos(c2) + B.* sin(c2)) ...
         -  exp(x1 * B) .* (B .* cos(c1) + B.* sin(c1)));

chi(1,1, 1) = x2 - x1; %Correct for the case j1 = j2 = 0!
end
