function phi_hest = chfun_heston( r, kappa, theta, V0, nu, rho, T, w, x)

%{
 This code computes the Characteristic Function for the Heston Model
%  Notation follows Fang and Oosterlee (2008), eq. 32, p. 8
% 
 Authors : Baldi Lanfranchi, Federico
         : La Cour, Peter

 Inputs 
            : r             - log price drift rate
            : kappa         - speed of mean reversion
            : theta         - mean (long run) volatility
            : V0            - initial volatility
            : nu            - volatility of the volatility (vol of vol)
            : rho           - correlation between Wiener processes (W1 and W2)
            : T             - time to maturity
            : w             - Vector of N evaluation intervals
 Optional   : x             - log asset price if not supplied, the last
                              term disappears
 Outputs    : phi_hest      - characteristic function values [0:N-1] vector
%}
%Default: use all generated paths!
if (~ exist ( 'x' , 'var' )) x = 0; end 

% D and G parameters: Fang and Oosterlee (2008) p. 8
D = sqrt( (kappa - 1i * rho * nu * w) .^ 2 + (w .^ 2 + 1i * w) * nu ^ 2 );
G = (kappa - 1i * rho * nu * w - D) ./ (kappa - 1i * rho * nu * w + D);

% Characteristic function for the Heston Model: 
% Fang and Oosterlee (2008) p. 8
phi_hest = exp( 1i * w * r * T + ...
           (V0 / nu ^ 2) * (1 - exp(-D * T)) ./ (1 - G .* exp(-D * T)) .* (kappa - 1i * rho * nu * w - D) + ...
           (kappa * theta / nu ^ 2 * (T * (kappa - 1i * rho * nu * w - D) - 2 * log( (1 - G .* exp(-D * T)) ./ (1 - G) ))) + ...
           1i * w * x); %note: if x = 0, this term disappears
      
  