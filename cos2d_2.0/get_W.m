function W = get_W(kappa, rho, eta, theta, r, dt, alpha, t, a1, b1, a2, b2, N, V)
  
    %to ease up notation
    pbma1 = pi / (b1 - a1); 
    pbma2 = pi / (b2 - a2);

    W = zeros(N); %empty matrice
    for j1 = 0:(N-1)
        for j2 = 0:(N-1)
            phi_A_ = phi_A(j1 * pbma1,  alpha * j2 * pbma2, kappa, rho, eta, theta, r, dt);
            %Note store at +1!!!
            W(j1 + 1, j2 + 1) = 0.5 * exp( - r * dt) * phi_A_ * V(j1 + 1, j2 + 1, t + 1);
        end
    end
end           
