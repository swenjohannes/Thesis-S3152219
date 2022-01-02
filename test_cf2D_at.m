parameter_set2_at

npoint = 32;

%H = 0.5;

%Truncation ranges
[c1, c2, ~]  = heston_cumulants_v1(r, q, kappa, theta, v0, eta, rho, T);
[a1, b1] = cos_truncation_range_v2(c1,c2,0,12); %Obtain a and b from cumulants
[a2, b2] = a2_b2(eta, theta, kappa, T, v0, 1e-4);
u1 = pi/(b1-a1)*(0:(npoint-1))'*2;
u2 = pi/(b2-a2)*(0:(npoint-1))*2;

%% Heston

fprintf('Standard Heston analytical: \t')
tic
[phi,A0hs,B1hs,B2hs] = cf2dHeston(u1,u2, x0(1),v0, r-q,kappa,theta,eta,rho,T);
toc()
figure(1)
surf(abs(phi))
%surf(real(phi))
%surf(imag(phi))
%surf(imag(B2hs))
title('Standard Heston analytical'), xlabel('u_2'), ylabel('u_1'), zlabel('|\phi(u)|')
view([100 20])

fprintf('Standard Heston Monte Carlo: \t')
tic
%phi_mc = cf2dHeston_mc(u1,u2, x0(1),v0, r-q,kappa,theta,eta,rho,T);
toc()
figure(2)
surf(abs(phi_mc))
%surf(real(phi_mc))
%surf(imag(phi_mc))
xlabel('u_2'), ylabel('u_1'), zlabel('|\phi(u)|')
view([100 20])

%% Rough Heston
% Volterra intergal equation
fprintf('Rough Heston Volterra: \t\t')
tic
[phi_hr,A0hr,B1hr,B2hr] = cf2dRoughHeston(u1,u2, x0(1),v0, r-q,kappa,theta,eta,rho,H,T,300);
toc()
figure(3)
surf(abs(phi_hr))
%surf(real(phi_hr))
%surf(imag(phi_hr))
%surf(imag(B2hr + 1i*ones(size(u1))*u2))
title(['Rough Heston Volterra, H=' num2str(H)]), xlabel('u_2'), ylabel('u_1'), zlabel('|\phi(u)|')
view([100 20])

% Monte Carlo
fprintf('Rough Heston Monte Carlo: \t')
tic
phi_hrmc = cf2dRoughHeston_mc(u1,u2, x0(1),v0, r-q,kappa,theta,eta,rho,H,T);
toc()
figure(4)
surf(abs(phi_hrmc))
%surf(real(phi_hrmc))
%surf(imag(phi_hrmc))
title(['Rough Heston Monte Carlo, H=' num2str(H)]), xlabel('u_2'), ylabel('u_1'), zlabel('|\phi(u)|')
view([100 20])

%% 1D
[p1d,B21d] = cfRoughHeston2(u1,r,v0,kappa,theta,eta,rho,H,T,nstep,x0(1));

figure(6)
plot(u1,imag(B21d),u1,imag(B2hr(:,1)),u1,imag(B2hs(:,1)))
