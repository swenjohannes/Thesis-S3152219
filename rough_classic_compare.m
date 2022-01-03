parameter_set2 %load set 3 to compare with ERKAN!

T = 1;
K = 80;
[S_c, V_c, time_c] = heston_mc(S0, v0, rho, kappa, theta, T, r, q, eta, npath, steps);
vanilla_prices(S_c, K, 1, r, T)

npath = 10;
steps= 1200;
[S_r, V_r, time_r] = rough_heston_mc4(S0, v0, rho, kappa, theta, T, r, q, eta, H, npath, steps);
vanilla_prices(S_r, 100, 1, r, T)
plot_SV(S_r, V_r, "rh", 1:1e4, 1:1e4)
cos2d("h", S0, K, 1000, v0, r, q, eta, theta, rho, kappa, 1, 100, 1) 
cos2d("rh", S0, K, 1000, v0, r, q, eta, theta, rho, kappa, 1, 100, 1, H, 0, 1.4423, 400)


cos2d("h", S0, K, 1000, v0, r, q, eta, theta, rho, kappa, 0.5, 100, 1) 

%cos2d("rh", S0, K, 1000, v0, r, q, eta, theta, rho, kappa, 1, 100, 1, 1, a2, b2, 500)

[c1, c2, ~]  = heston_cumulants_v1(r, q, kappa, theta, v0, eta, rho, T);
[a1, b1] = cos_truncation_range_v2(c1,c2,0,12); %Obtain a and b from cumulants
[a2, b2] = a2_b2(eta, theta, kappa, T, v0, 1e-4);
Nobs = 1;
dt = T/Nobs;
k = 0:(N-1);        %Pre compute A and B2 parts of the characteristic equation!

[phi_Ap, B2p]  = phi_A_B2_heston(k, k, kappa, rho, eta, theta, r, q, dt, a1, b1, a2, b2);
[phi_Am, B2m]  = phi_A_B2_heston(k, -k, kappa, rho, eta, theta, r, q, dt, a1, b1, a2, b2);

[phi_Ap_r, B2p_r]  = phi_A_B2_rough_heston(k, k, kappa, rho, eta, theta, r, q, a1, b1, a2, b2, 1, dt, 500);
[phi_Am_r, B2m_r]  = phi_A_B2_rough_heston(k, -k, kappa, rho, eta, theta, r, q, a1, b1, a2, b2, 1, dt, 500);
