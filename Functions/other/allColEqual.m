function test = allColEqual(A, tol)

%Set default tolerance to 1e-6
if (~ exist ( 'tol' , 'var' )) tol = 1e-6; end

test = abs(diff(A, 1, 2)) < tol;

if sum(~test, 'all') == 0
    disp("All columns are equal")
else 
    disp("At least one or more columns is unequal")
end