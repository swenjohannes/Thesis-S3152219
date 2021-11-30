function H = get_H(kappa, rho, eta, dt, alpha, a, b, N)

%{
    Description: Computes a matrix containing all H(k2, j2, j1) values. 
    
    Parameters:
      kappa:    [1x1 real] Heston parameter
      rho:      [1x1 real] Heston parameter
      eta:      [1x1 real] Heston parameter
      dt:       [1x1 real] time between observations
      alpha:    [1x1 integer] either 1 or -1 to compute positive or
                negative part
      a:        [1x1 real] first cosine argument
      b:        [1x1 real] second cosine argument

    Output: 
      H:        [N x N x N] matrix indexed by (k2, j2, j1)
    
    References:
      - Moir 4.31 
%}

    bma = b - a;  
    pbma = pi / bma;  
    
    k2 = 0:(N-1);           %Compute K2 values in vector vorm
    H(:,:, N) = zeros(N);   %empty matrice
    for j1 = 0:(N-1)
        for j2 = 0:(N-1)                        %Compute B2                        
                B = B2(j1 * pbma, alpha * j2 * pbma, kappa, rho, eta, dt);
                %Compute Moir 4.32, chi_H returns a vector over all k2's
                H(:, j2 + 1, j1 + 1) =  2 / bma * exp(alpha * 1i * j2 * -a *pbma) ...
                                     .* chi_H(a, b, B, 0:(N-1), a, b);
        end
    end
end           
