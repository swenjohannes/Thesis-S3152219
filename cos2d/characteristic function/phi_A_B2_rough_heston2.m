function [phi_A, B2] = phi_A_B2_rough_heston2(u1, u2, kappa, rho, eta, theta, r, q, a1, b1, a2, b2, alpha, T, N)
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

%u1 = 1:10; u2 =1:10; alpha = 1; N = 160;

pbma1 = pi / (b1 - a1);
pbma2 = pi / (b2 - a2);
u1 = u1' * pbma1; %make transpose!
dt = T / N;
M = length(u1); 

% To define the Volterra integral equation:
c1 = - 0.5 *(u1 .^2 + 1i * u1);
beta = kappa - 1i * rho * eta * u1; 
c3 = 0.5*eta ^ 2;

[B2, phi_A] = deal(zeros([M, M]));
for j = u2
    j = 5;
    g = @(t) 1i * j * pbma2 * t.^-alpha / gamma(1-alpha);
    f = @(w) (c1 - beta .*w + c3*w.^2);
    [psi, Dalpha_psi, Dalpha_psi2] = SolveVIE2(f, g,  1i * j * pbma2, alpha,T, N, M); %solve VIE
%     close all;
    plot_imag(Dalpha_psi)
    plot_imag(Dalpha_psi2)
    % Integrate to get the characteristic function
%       B2(:, abs(j) + 1) = sum((Dalpha_psi(:, 1:end-1) + Dalpha_psi(:,2:end))/2, 2) *dt  ...
%                              + 1i * j * pbma2; %Don't know why this needs to be added.

    Dalpha_psi2(:, 1)
    B2(:, abs(j) + 1) = sum((Dalpha_psi2(:, 1:end-1) + Dalpha_psi2(:,2:end))/2, 2) *dt ...
                         + 1i * j * pbma2; %Don't know why this needs to be added.

    phi_A(:, abs(j) + 1) = exp(1i * u1 * (r - q) * T + ...
                    kappa * theta * sum((psi(:, 1:end-1) + psi(:, 2:end))/2, 2) *dt);
end
end
% 
% alpha = 0.6
% 0.9 ^(1 - alpha) 
% 1 - (0.5 / 12) ^-alpha / gamma(1 - alpha)