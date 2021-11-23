function [C_carr_fft, ku] = cm_vannila_fft(psi, a, N, eta)
%using the approach of Carr-Madan 1999 (FFT)

%Default parameters
if (~ exist ( 'a' , 'var' )) a = 1.2; end 
if (~ exist ( 'N' , 'var' )) N = 2^10; end 
if (~ exist ( 'eta' , 'var' )) eta = 0.5; end 

%Initialization
lambda = (2 * pi) / ( N * eta);
b =  N * lambda / 2;
u = 1:N;
ku = -b + lambda * (u - 1);
vj = eta * ( u - 1);
w = exp(1i * b * vj) .*  psi(vj);

%Calculate G term
g1 = g(vj(1), ku(1), psi);
gN = g(vj(N), ku(N), psi); 
gterm = 0.5 * (g1 + gN);

%Finaly compute prices
C_carr_fft = exp(-a * ku) ./ pi .* real(eta * (fft(w) - gterm));
end