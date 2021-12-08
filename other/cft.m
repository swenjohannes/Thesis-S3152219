function F = cft(X,n)
% COS Fourier transformation

m = size(X,1);
Y = linspace(0,1,m);
dY = Y(2) - Y(1);
K = pi*(0:n-1)';

F = cos(K * Y) * X * (2 * dY);