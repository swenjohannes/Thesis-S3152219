function vec = vectorize_fun(v, psi_N_)
n = length(v);
vec = ones(1, length(v)); %create empty vector
for j = 1:n
    vec(j) = psi_N_(j);
end
end