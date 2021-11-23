function A = get_A(Hp, Hm, Wp, Wm, N)

A = zeros(N); %Empty frame

for j1 = 1:N
    A(j1, :) = Hp(:, :, j1) * Wp(j1, :)' + Hm(:, :, j1) * Wm(j1, :)' ;
end
end