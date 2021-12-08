function X = icft2d(F,m,n)
% Inverse 2D COS Fourier transformation

X1 = icft(F',n)';
X = icft(X1,m);
