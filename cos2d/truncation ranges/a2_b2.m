function [a2, b2] = a2_b2(eta, theta, kappa, T, v0, tol)

%{
    Description: Computes a matrix containing all H(k2, j2, j1) values. 
    
    Parameters:
      eta:      [1x1 real] Heston parameter
      theta:    [1x1 real] Heston parameter
      kappa:    [1x1 real] Heston parameter
      T:        [1x1 real] Time to maturity
      v0:       [1x1 real] Initial variance
      tol:      [1x1 real] Tolerance level
    
    Output: 
      a2:       [1x1 real] lower truncation bound of variance
      b2:       [1x1 real] upper truncation bound of variance    
    References:

%}

C0 = eta ^ 2 * (1 - exp(- kappa * T)) / ( 4 * kappa); %Compute the constant
d = 4 * kappa * theta / eta ^ 2;                      %Degrees of freedom
lambda = 4 * kappa * exp(- kappa * T) * v0 ...        %Non-centrality ..
        / (eta ^ 2 * (1 - exp(- kappa * T)));         %parameter

a2 = ncx2inv(tol, lambda,d) * C0;
b2 = ncx2inv(1 - tol, lambda,d) * C0;
end