function F = get_F2(x0, v0, phi_Ap, phi_Am, B2p, B2m, a1, b1, a2, b2, N)

%{
    Description: Computes the F = 1/2 (Fp + Fm)
    
    Parameters:  
                x0         [1x1 real] initial log spot price
                v0         [1x1 real] initial variance
                Ap, Am:    [j1 x j2 complex] Positive/negative A part of phi
                B2p, B2m:  [j1 x j2 complex] Positive/negative B2 part of phi
                a1, b1:    [1x1 real] Cosine arguments of log spot price
                a2, b2:    [1x1 real] Cosine arguments of volatility
                N          [1x1 real] Truncation argument

    Output: 
      A:        [j1 x k2] matrix containing weighted_sum(Hp * Wp + Hm * Wm)
    
    References:

%}
phi_p = phi(x0, v0, phi_Ap, B2p, a1, b1, N); %Obtain positive charfn
phi_m = phi(x0, v0, phi_Am, B2m, a1, b1, N); %Obtain negative charfn

k1 = (0:(N-1))' * pi * a1 / (b1 - a1); %Vector compute 
k2 = (0:(N-1)) * pi * a2 / (b2 - a2);  %Vector compute 

F = 1/2 * real(phi_p .* exp(-1i * (k1 + k2))  ...
             + phi_m .* exp(-1i * (k1 - k2))); 
end
