function phi_hest = heston_char_fn( r, kappa, theta, V0, eta, rho, T, omega, x)

%{
 This code computes the Characteristic Function for the Heston Model
 Notation follows Fang and Oosterlee (2008), eq. 32, p. 8

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
            : a             - Cosine argument (lower truncation bound) 
            : b             - Cosine argument (upper truncation bound)
            : k             - Vector of N evaluation intervals
 Optional   : x             - log asset price if not supplied, the last
                              term disappears
 Outputs    : phi_hest      - characteristic function values [0:N-1] vector
%}

if (~ exist ( 'x' , 'var' )) x = 0; end 

% D and G parameters: Fang and Oosterlee (2008) p. 8
D = sqrt( (kappa - 1i * rho * eta * omega) .^ 2 + (omega .^ 2 + 1i * omega) * eta ^ 2 );
G = (kappa - 1i * rho * eta * omega - D) ./ (kappa - 1i * rho * eta * omega + D);

% Characteristic function for the Heston Model: 
% Fang and Oosterlee (2008) p. 8
phi_hest = exp( 1i * omega * r * T + ...
           (V0 / eta ^ 2) * (1 - exp(-D * T)) ./ (1 - G .* exp(-D * T)) .* (kappa - 1i * rho * eta * omega - D) + ...
           (kappa * theta / eta ^ 2 * (T * (kappa - 1i * rho * eta * omega - D) - 2 * log( (1 - G .* exp(-D * T)) ./ (1 - G) ))) + ...
           i * omega * x); %note: if x = 0, this term disappears
      
  