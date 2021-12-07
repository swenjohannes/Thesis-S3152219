function chi = chi_H(x1, x2, B2, a, b, N)
%{
    Description: Computes the integral of H analytical
    
    Parameters:
      x1:       [1x1 real] Lower integration range
      x2:       [1x1 real] Upper integration range
      B2:       [j1 x j2]  Matrix of B2 part of char eq
      a:        [1x1 real] Lower cosine argument
      b:        [1x1 real] Upper cosine argument
      N:        [1x1 real] Truncation argument

    Output: 
      chi:      [j1 x j2 x k2] matrix containing all chi_H values
    
    References:

%}

k = zeros(1); %Create 1x1 element
k(1, 1, N) = 0; %Add 3D dimension
k(1, 1, :) =  (0:(N-1)) * pi / (b - a); %Give proper values

%Arguments of cos/sin
c1 = k * (x1 - a);
c2 = k * (x2 - a);

chi = 1 ./ (B2 .^2 + k .^ 2) ...
    .* ( exp(x2 * B2) .* (B2 .* cos(c2) + B2.* sin(c2)) ...
         -  exp(x1 * B2) .* (B2 .* cos(c1) + B2.* sin(c1)));

chi(1,1, 1) = x2 - x1; %Correct for the case j1 = j2 = 0!
end
