function [phi_A, B2] = phi_A_B2_rough_heston(u1, u2, kappa, rho, eta, theta, r, q, a1, b1, a2, b2, alpha, T, N)
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
  q:        [1x1 real] Heston parameter
  tau:      [1x1 real] Time delta between observations
  a1, b1:   [1x1 real] Cosine arguments of log spot price
  a2, b2:   [1x1 real] Cosine arguments of variance
  N         [1x1 real] Truncation argument

Output: 
  phi_A:    [N x M] matrix containing all phi_A values for the input w1,
             w2. Note the function does wj * pi/(bj -aj)!
  B2:       [N x M] matrix containing all B2 values for the input w1,
             w2. Note the function does wj * pi/(bj -aj)!
References:

%}

pbma1 = pi / (b1 - a1);
pbma2 = pi / (b2 - a2);
u1 = u1' * pbma1; %make transpose!
dt = T / N;
M = length(u1); 

% Define the Volterra integral equation:
c1 = - 0.5 *(u1 .^2 + 1i * u1);
beta = kappa - 1i * rho * eta * u1; 
c3 = 0.5*eta ^ 2;
f = @(w) (c1 - beta .*w + c3*w.^2);

[B2, phi_A] = deal(zeros([M, M]));
for j = u2
    [psi,Dalpha_psi] = SolveVIE(f, 1i * j * pbma2, alpha,T, N, M); %solve VIE
 
    % Integrate to get the characteristic function
    B2(:, abs(j) + 1) = sum((Dalpha_psi(:, 1:end-1) + Dalpha_psi(:,2:end))/2, 2) *dt ...
                        + 1i * j * pbma2; %Don't know why this needs to be added.
    phi_A(:, abs(j) + 1) = exp(1i * u1 * (r - q) * T + ...
                    kappa * theta * sum((psi(:, 1:end-1) + psi(:, 2:end))/2, 2) *dt);
end
end

