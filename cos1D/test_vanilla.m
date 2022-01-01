% Test param
parameter_set2_at

pc = 1;
nstep = 256;
npath = 1e4;
K = [   1.094900000000000
   1.099200000000000
   1.103400000000000
   1.107700000000000
   1.111900000000000
   1.116200000000000
   1.120400000000000
   1.124600000000000
   1.128900000000000
   1.133100000000000
   1.137400000000000
   1.141600000000000
   1.145900000000000
   1.150100000000000
   1.154400000000000
   1.158600000000000
   1.162800000000000]';
%

%% Heston
disp('Heston')
tic
for i = 1:1
Vcos = hsVanilla_cos(K,pc,T, S0,r,q,v0,kappa,theta,eta,rho);
end
Vcos
toc()
tic
Vmc = hsVanilla_mc(K,pc,T, S0,r,q,v0,kappa,theta,eta,rho, npath,nstep)
toc()

%% Rough Heston
disp('Rough Heston')

tic
for i = 1:1
    Vcos = hrVanilla_cos(K,pc,T, S0,r,q,v0,kappa,theta,eta,rho,H)
end
toc()

tic
Vmc = hrVanilla_mc(K,pc,T, S0,r,q,v0,kappa,theta,eta,rho,H, npath,nstep)
toc()

