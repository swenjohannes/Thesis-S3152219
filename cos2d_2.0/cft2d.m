function F = cft2d(X,m,n)
% 2D COS Fourier transformation

F1 = cft(X,m);
F = cft(F1',n)';