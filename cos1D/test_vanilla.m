% Test param
parameter_set2_at
pc = 1;
nstep = 256;
npath = 1e4;
K = (50:10:150);
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
    Vcos = hrVanilla_cos(K,pc,T, S0,r,q,v0,kappa,theta,eta,rho,H, 128)
end
toc()

tic
Vmc = hrVanilla_mc(K,pc,T, S0,r,q,v0,kappa,theta,eta,rho,H, npath,nstep)
toc()

