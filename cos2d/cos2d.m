function V0 = cos2d(model, S0, K, L, v0, r, q, eta, theta, rho, kappa, T, N, Nobs, alpha)

%{
    Description: Computes the barrier price under the Heston model, using
                 the two dimensional COS method.
    
    Parameters:
      model     [string]   Model argument
      S0:       [1x1 real] Stock price
      K:        [1x1 real] Strike price
      L:        [1x1 real] Barrier level 
      v0:       [1x1 real] Initial variance
      r:        [1x1 real] Riskfree interest rate
      q:        [1x1 real] Divident yield
      eta:      [1x1 real] Volatility of volatility
      theta:    [1x1 real] Longrun volatility level
      rho:      [1x1 real] Correlation coefficient
      kappa:    [1x1 real] Mean reversion speed volatility 
      T:        [1x1 real] Time to maturity
      N:        [1x1 real] Number of cosine arguments
      N_obs:    [1x1 real] Number of observations, can be multiple!
    
Output: 
      V0:       [1x1 real] Barrier price
    References:

%}
%% Initialization
x0 = log(S0 / K);                           %Log strike-spot
h = log(L / K);                             
dt = T/ Nobs;                               %Time steps

%% Truncation ranges
[c1, c2, ~]  = heston_cumulants_v1(r, q, kappa, theta, v0, eta, rho, T);
[a1, b1] = cos_truncation_range_v2(c1,c2,0,12); %Obtain a and b from cumulants
[a2, b2] = a2_b2(eta, theta, kappa, T, v0, 1e-4);

%% Pre compute A and B2 parts of the characteristic equation!
k = 0:(N-1);        
if model == "h"
    alpha = NaN; %Optional parameter
    [phi_Ap, B2p]  = phi_A_B2_heston(k, k, kappa, rho, eta, theta, r, q, dt, a1, b1, a2, b2);
    [phi_Am, B2m]  = phi_A_B2_heston(k, -k, kappa, rho, eta, theta, r, q, dt, a1, b1, a2, b2);
elseif model == "rh"
    [phi_Ap, B2p]  = phi_A_B2_rough_heston(k, k, kappa, rho, eta, theta, r, q, a1, b1, a2, b2, 1, dt, 500);
    [phi_Am, B2m]  = phi_A_B2_rough_heston(k, -k, kappa, rho, eta, theta, r, q, a1, b1, a2, b2, alpha, dt, 500);
else error('Please supply a correct model. Options are "h" or "rh". ')   %Succesful!
end

%% COS 2D
V(:, :, Nobs) = zeros(N);           %Create empty V dataframe
V(:, :, Nobs) = V_T(h, a1, b1, a2, b2, K, N); %Store V(T_M) at M
if Nobs > 1
    Hp = get_H(1, B2p, a2, b2, N);                         %positive part of H
    Hm = get_H(-1, B2m, a2, b2, N);                        %negative part of M
    Mp = get_M(a1, h, a1, b1, N);                           %Obtain M

    for t = (Nobs - 1):-1:1
        Wp = 0.5 * exp(-r * dt) * phi_Ap .* V(:, :, t + 1); %W positive
        Wm = 0.5 * exp(-r * dt) * phi_Am .* V(:, :, t + 1); %W minus

        Wp(:,1) = 0.5 * Wp(:,1);  % Set first j2 element of both W matrices to
        Wm(:,1) = 0.5 * Wm(:,1);  % a half to compute the weighted sum!

        A = reshape(sum(Hp .* Wp + Hm .* Wm, 2), [N N]); %Take rowsums and remove 3D!
        A(1, :) = A(1, :) * 0.5; %Weight first row by a half!
        V(:, :, t) = real(Mp * A); %Compute the weighted sum in matrix form
    end
end

%% Final step: compute V0
F = get_F(x0, v0, phi_Ap, phi_Am, B2p, B2m, a1, b1, a2, b2, N);

F(1, :) = F(1, :) * 0.5; %weight first row a by a half
F(:, 1) = F(:, 1) * 0.5; %weigth first col by a half
V0 = exp(-r * dt)  * sum(F .* V(:, : , 1), 'all'); %Calc V0 as 2D weigthed sum
end
