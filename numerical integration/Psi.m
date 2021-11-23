function res = Psi(v, ch_fun, r, T, a)
    res = exp(-r * T) .* ch_fun(v - (a + 1) * 1i)  ... 
            ./ (a ^ 2 + a - v .^ 2 + 1i * (2 * a + 1) * v);
end