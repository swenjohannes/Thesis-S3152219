function X = icft(F,n)
% Inverse COS Fourier transformation

m = size(F,1);
Y = linspace(0,1,n)';
K = pi*(0:m-1);

F(1,:) = F(1,:)/2;
X = cos(Y * K) * F;

