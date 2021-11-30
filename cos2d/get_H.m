function H = get_H(kappa, rho, eta, dt, alpha, a, b, N)
    % Calculate H for all k2, j2, j1 using formula Moir 4.31

    bma = b - a;  
    pbma = pi / bma;  
    
    k2 = 0:(N-1);
    H(:,:, N) = zeros(N); %empty matrice
    for j1 = 0:(N-1)
        for j2 = 0:(N-1)
                B = B2(j1 * pbma, alpha * j2 * pbma, kappa, rho, eta, dt);
                H(:, j2 + 1, j1 + 1) =  2 / bma * exp(alpha * 1i * j2 * -a *pbma) ...
                                     .* chi_H(a, b, B, k2, a, b);
        end
    end
end           

% 
% Check for equality of chi_H and integral
% 
% bma = b - a;  
% pbma = pi / bma;  
% alpha = 1;
% j1 = 20;
% j2 = 1;
% k2 = 1;
%  
% B = B2(j1 * pbma, alpha * j2 * pbma, kappa, rho, eta, dt)
% chi_H(a, b, B, k2, a, b)
% 2 / bma * exp(alpha * 1i * j2 * -a *pbma) ...
%                      .* chi_H(a, b, B, k2, a, b)
% 
%  chi_H(a, b, B, k2, a, b)
% integrand = @(y) exp(y * B) .* cos(k2 * pbma * (y - a));
% integral(integrand, a, b)