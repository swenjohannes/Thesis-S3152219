function phi = phi(x, v, phi_A, B2, a1, b1, N)
%{
Description: Calculates the A part of the heston characteristic function.

Parameters:
  x:        [1x1 real] Initial log stock price
  v:        [1x1 real] Initial variance
  phi_A:    [N x N]    Matrix of pre-computed phi_A values
  B2:       [N x N]    Matrix of pre-computed B2 values
  a1, b1:   [1x1 real] Cosine arguments of log spot price
  N         [1x1 real] Truncation argument

Output: 
  phi:        [N x N] matrix containing all phi values

References:

%}

w1 = (0:(N-1))' * pi / (b1 - a1);  %Turn w1 in w1* and make a column vector
phi = exp(1i * w1 * x + B2 * v) .* phi_A; %Compute phi
end