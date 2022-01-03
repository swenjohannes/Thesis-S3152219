function iv = get_iv_hr(K,pc,T, S0,r,q,x)
disp(x)
v0 = x(1);
kappa = x(2);
theta = x(3);
eta = x(4);
rho = x(5);
H = x(6);

price  = zeros(length(T),length(K));
iv  = zeros(size(price));

for i = 1:length(T)
    price(i, :) = hrVanilla_cos(K,pc,T(i), S0,r,q,v0,kappa,theta,eta,rho, H);
    iv(i,:) = bsimpvol(S0,K,r,q,T(i),pc,price(i,:));
end

end