function C_carr = cm_vannila(psi, k, a)
%Numerical integration
integrand = @(v) exp(-1i * v * k) .* psi(v);
C_carr = exp(-a * k) / pi * real(integral(integrand, 0, 100));
end