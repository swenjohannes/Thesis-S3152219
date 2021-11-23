function M = mkj_mat(N, x1, x2, a, b)
M = zeros(N);
mkj_ = @(k, j) mkj(k, j, x1, x2, a, b);
for k = 1:N
    for j = 1:N
        M(k, j) = mkj_(k - 1, j - 1); %k and j start from 0, but are stored at loc 1
    end
end
