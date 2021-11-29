function Sum = weighted_sum(x, y)
if (~ exist ( 'y' , 'var' )) y = 1; end  %To use the function 1D 
    f = x .* y;
    f(1) = 0.5 * f(1);
    Sum = sum(f);
end