function J = get_H(kappa, rho, eta, dt, alpha, a, b, N)
    %to ease up notation
    bma = b - a;
    pbma = pi / bma;  
    B2_ = @(w1, w2) B2(w1 * pbma, w2 * pbma, kappa, rho, eta, dt);

    J(:,:, N) = zeros(N); %empty matrice

    k2 = 0:(N-1);
    for j1 = 0:(N-1)
        for j2 = 0:(N-1)
                J(:, j2 + 1, j1 + 1) =  2 / (b - a) * exp(alpha * 1i * j2 * -a *pbma) ...
                    .* chi_H(a, b, B2_(j1, alpha * j2), k2, a, b);
        end
    end
end           
