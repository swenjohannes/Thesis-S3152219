function phi = chfun_heston(v0, kappa, rho, eta, theta, r, q, tau, w)
%{
Description: Calculates the A part of the heston characteristic function.

Parameters:
  w:       [N x 1] vector of w values
  kappa:    [1x1 real] Heston parameter
  rho:      [1x1 real] Heston parameter
  eta:      [1x1 real] Heston parameter
  theta:    [1x1 real] Heston parameter
  r:        [1x1 real] Heston parameter
  q:        [1x1 real] Heston parameter
  tau:      [1x1 real] Time delta between observations

Output: 
  phi:    matrix containing char func results
References:

%} 

eta2 = eta ^ 2; %Notation!
beta = kappa - 1i * rho * eta * w; 
D = sqrt(beta .^ 2 + eta2 * w .* (1i + w));

h = (beta - D) ./ (beta + D);
hemdt = h  .* exp (-D * tau); %To ease up calculations
B2 = (beta - D - (beta + D) .* hemdt) ./ (eta2 * (1 - hemdt)); %Get B2
A = 1i * w * (r - q) * tau + ...  %Compute phi A
    kappa * theta / eta2 * ((beta - D) * tau - 2 * log((hemdt - 1) ./ (h - 1)));
phi = exp(B2 * v0 + A);

end