function price = hsVanilla_cos(K,pc,T, S0,r,q,v0,kappa,theta,eta,rho)

% % Test param
% parameter_set2
% pc = -1;
% %

S0 = S0(:)';
K = K(:)';
x = log(S0 ./ K);
N = 256;

[c1, c2, ~]  = heston_cumulants_v1(r, q, kappa, theta, v0, eta, rho, T);
[a, b] = cos_truncation_range_v2(c1,c2,0,12); %Obtain a and b from cumulants
w = pi/(b-a)*(0:(N-1))';

%Obtain trucation ranges from cumulants
%chfun_heston( r, kappa, theta, V0, nu, rho, T, w, x)
phi = chfun_heston(r-q, kappa, theta, v0, eta, rho, T, w);
phi(1, :) = 0.5 * phi(1, :); %weight first term by a half

%Uk = 2 / (b - a) *(chi_coef(k, 0, b, a, b) - psi_coef(k, 0, b, a, b));
if pc >= 0 % call
    Uk = 2 / (b - a) *(chi(w, 0, b, a, b) - psi(w, 0, b, a, b));
else % put
    Uk = 2 / (b - a) *(-chi(w, a, 0, a, b) + psi(w, a, 0, a, b));
end

% F = (phi .* Uk) .* exp(1i * w * (x - a));
% price = K * exp(-r * T) .* real(sum(F, 1));

F = (phi .* Uk).' * exp(1i * w * (x - a));
price = K * exp(-r * T) .* real(F);

price(price<0) = 0;
end
function xout = chi(w,c,d,a,b)
    c1 = ( cos(w*(d-a)) + w .* sin(w*(d-a)) ) * exp(d);
    c2 = ( cos(w*(c-a)) + w .* sin(w*(c-a)) ) * exp(c);
    xout = (c1 - c2)./(1+w.^2);
end
function pout = psi(w,c,d,a,b)
    pout = (d-c) * ones(size(w));
    i = ( w~= 0 );
    s1 = sin(w(i)*(d-a));
    s2 = sin(w(i)*(c-a));
    pout(i) = (s1 - s2)./w(i);
end