function J = get_H(kappa, rho, eta, dt, alpha, a, b, N)
    bma = b - a;  %to ease up notation
    pbma = pi / bma;  
    
    J(:,:, N) = zeros(N); %empty matrice

    k2 = 0:(N-1);
    for j1 = 0:(N-1)
        for j2 = 0:(N-1)
                B = B2(j1 * pbma, alpha * j2 * pbma, kappa, rho, eta, dt);
                J(:, j2 + 1, j1 + 1) =  2 / bma * exp(alpha * 1i * j2 * -a *pbma) ...
                                     .* chi_H(a, b, B, k2, a, b);
        end
    end
end           


% 
% Check for equality of chi_H and integral
% 
% alpha = 1;
% j1 = 7;
% j2 = 1;
% k2 = 1;
%  
% B = B2(j1 * pbma, alpha * j2 * pbma, kappa, rho, eta, dt);
%   2 / bma * exp(alpha * 1i * j2 * -a *pbma) ...
%                      .* chi_H(a, b, B, k2, a, b)

%  chi_H(a, b, B, k2, a, b)
% integrand = @(y) exp(y * B) .* cos(k2 * pbma * (y - a));
% integral(integrand, a, b)