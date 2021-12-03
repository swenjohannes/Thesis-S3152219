function phi_A_mat = get_phi_A(kappa, rho, eta, theta, r, dt, a1, b1, a2, b2, N, alpha)
%to ease up notation
pbma1 = pi / (b1 - a1);  
pbma2 = pi / (b2 - a2); 

phi_A_mat = zeros(N);
for j1 = 0:(N-1)
   for j2 = 0:(N-1)
       phi_A_mat(j1 + 1, j2 + 1) = phi_A(j1 * pbma1,  alpha * j2 * pbma2, kappa, rho, eta, theta, r, dt);
   end
end
end