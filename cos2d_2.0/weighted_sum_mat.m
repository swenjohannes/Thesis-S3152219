function Sum = weighted_sum_mat(A)
N = length(A);
Sum = 0.5 * weighted_sum(A(1, :)); %first times 1/2!
for i = 2:N
    Sum = Sum + weighted_sum(A(i, :));
end

%Example usage: 245.25 is correct!
%U = [1,2,3; 4,5,6; 7,8,9]
%weighted_sum_mat(U,U)