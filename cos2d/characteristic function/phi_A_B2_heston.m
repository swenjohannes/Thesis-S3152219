function [phi_A, B2] = phi_A_B2_heston(w1, w2, kappa, rho, eta, theta, r, q, tau)
%{
Description: Calculates the A part of the heston characteristic function.

Parameters:
  w1:       [N x 1] vector of w1 values.
  w2:       [M x 1] vector of w1 values.
  kappa:    [1x1 real] Heston parameter
  rho:      [1x1 real] Heston parameter
  eta:      [1x1 real] Heston parameter
  theta:    [1x1 real] Heston parameter
  r:        [1x1 real] Heston parameter
  q:        [1x1 real] Heston parameter
  tau:      [1x1 real] Time delta between observations
  N         [1x1 real] Truncation argument

Output: 
  phi_A:    [N x M] matrix containing all phi_A values for the input w1,
             w2. Note the function does wj * pi/(bj -aj)!
  B2:       [N x M] matrix containing all B2 values for the input w1,
             w2. Note the function does wj * pi/(bj -aj)!
References:

%}

eta2 = eta ^ 2; %Notation!
beta = kappa - 1i * rho * eta * w1; 
D = sqrt(beta .^ 2 + eta2 * w1 .* (1i + w1));

h = (beta - D - 1i * w2 * eta2) ./ (beta + D - 1i * w2 * eta2);
hemdt = h  .* exp (-D * tau); %To ease up calculations

B2 = (beta - D - (beta + D) .* hemdt) ./ (eta2 * (1 - hemdt)); %Get B2

phi_A = exp(1i * w1 * (r - q) * tau + ...  %Compute phi A
    kappa * theta / eta2 * ((beta - D) * tau - 2 * log((hemdt - 1) ./ (h - 1))));

end