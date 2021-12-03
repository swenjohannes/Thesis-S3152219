function phi_A = phi_A(w1, w2, kappa, rho, eta, theta, r, tau, a1, b1, a2, b2)
%{
Description: Calculates the A part of the heston characteristic function.

Parameters:
  w1:       [N x 1] vector of w1 values. Should not be w1*!
  w2:       [M x 1] vector of w1 values. Should not be w1*!
  kappa:    [1x1 real] Heston parameter
  rho:      [1x1 real] Heston parameter
  eta:      [1x1 real] Heston parameter
  theta:    [1x1 real] Heston parameter
  r:        [1x1 real] Heston parameter
  tau:      [1x1 real] Time delta between observations
  a1, b1:   [1x1 real] Cosine arguments of log spot price
  a2, b2:   [1x1 real] Cosine arguments of variance

Output: 
  phi_A:    [N x M] matrix containing all phi_A values for the input w1,
             w2. Note the function does wj * pi/(bj -aj)!
References:

%}

w1 = w1' * pi / (b1 - a1);  %Turn w1 in w1* and make a column vector
w2 = w2  * pi / (b2 - a2);

eta2 = eta ^ 2; %Notation!
beta = kappa - 1i * rho * eta * w1; 
D = sqrt(beta .^ 2 + eta2 * w1 .* (1i + w1));

h = (beta - D - 1i * w2 * eta2) ./ (beta + D - 1i * w2 * eta2);
hemdt = h  .* exp (-D * tau); %To ease up calculations

phi_A = exp(1i * w1 * r * tau + ... %Compute phi A
    kappa * theta / eta2 * ((beta - D) * tau - 2 * log((hemdt - 1) ./ (h - 1))));
end