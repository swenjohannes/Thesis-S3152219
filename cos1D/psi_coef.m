function [psi_k] = psi_coef(k, x1, x2, a, b)

%{
 Inputs : a             - Cosine argument (lower truncation bound) 
        : b             - Cosine argument (upper truncation bound)
        : x1            - x1 argument
        : x2            - x2 argument
        : k             - FFT Nx1 vector of k [0:N-1]

Outputs : psi_k         - Nx1 vector of first psi cosine coefficients 
%}


bma = b - a;
pbma = pi / bma;
k_star = pbma * k;

%Arguments of cos/sin
c1 = k_star * (x1 - a);
c2 = k_star * (x2 - a);

%compute psi
psi_k = 1 ./ k_star .* (sin(c2) - sin(c1));
psi_k(1) = x2 - x1;
end