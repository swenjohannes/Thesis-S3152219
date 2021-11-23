function Sum = weighted_sum(x, y)
    f = x .* y;
    f(1) = 0.5 * f(1);
    Sum = sum(f);
end