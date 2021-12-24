% Test param
parameter_set2
pc = 1;
nstep = 200;
npath = 1e4;
K = 90;
%

%% Heston
disp('Heston')
tic
Vcos = hsVanilla_cos(K,pc,T, S0,r,q,v0,kappa,theta,eta,rho)
toc()
tic
Vmc = hsVanilla_mc(K,pc,T, S0,r,q,v0,kappa,theta,eta,rho, npath,nstep)
toc()

%% Rough Heston
disp('Rough Heston')

tic
Vmc = hrVanilla_mc(K,pc,T, S0,r,q,v0,kappa,theta,eta,rho,H, npath,nstep)
toc()
