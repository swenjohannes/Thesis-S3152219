function phi = chfun_rough_heston(r, q, kappa, theta, v0, eta, rho, T, w, H, N)
    %N_int = 160; w = k * pi / (b -a); H = 0.5; V0 = v0

    alpha = H + 0.5;
     % Define the Volterra integral equation:
    c1 = - 0.5 *(w .^2 + 1i * w);
    beta = kappa - 1i * rho * eta * w; 
    c3 = 0.5*eta ^ 2;
    f = @(y) (c1 - beta .*y + c3*y.^2);

    [psi,Dalpha_psi] = SolveVIE(f, 0, alpha,T, N, length(w)); %solve VIE SolveVIE(f, y0, alpha, T, N, M)

    dt = T / N;
    phi = exp( 1i * w * (r - q) * T ...
              + theta*kappa*sum((psi(:, 1:end-1) + psi(:,2:end))/2, 2) *dt...
               + v0.*sum((Dalpha_psi(:, 1:end-1) + Dalpha_psi(:,2:end))/2, 2) *dt);
end

