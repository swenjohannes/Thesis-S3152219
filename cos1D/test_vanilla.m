% Test param
parameter_set2
pc = 1;
nstep = 200;
npath = 1e4;
K = (80:10:100);
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
%Vmc = hrVanilla_mc(K,pc,T, S0,r,q,v0,kappa,theta,eta,rho,H, npath,nstep)
toc()

tic
%Vmc = hrVanilla_cos(K,pc,T, S0,r,q,v0,kappa,theta,eta,rho,H)
toc()
