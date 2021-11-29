function A = get_A(Hp, Hm, Wp, Wm, N)

A = zeros(N); %Empty frame

for j1 = 1:N
    for k2 = 1:N
        A(j1, k2) = weighted_sum(Hp(k2, :, j1) .* Wp(j1, :) + Hm(k2, :, j1) .* Wm(j1, :));
    end
end


% Test
% j1 = 6;
% k2 = 1;
% 
% p = Hp(k2, :, j1) .* Wp(j1, :);
% m = Hm(k2, :, j1) .* Wm(j1, :);
% weighted_sum(p,m)
% 

% Test 2 : weighted sum checken
% j1 = 10;
% k2 = 5;
% 
% phi_A_ = @(j1, j2) phi_A(j1 * pbma, j2 * pbma, kappa, rho, eta, theta, r, dt);
% 
% Sum = 0.5 * (0.5 * V(j1, 1, 2) .* (Hp(k2, 1, j1) .* phi_A_(j1 - 1, 0) + Hm(k2, 1, j1) .* phi_A_(j1 - 1, 0)) );
% for j2 = 2:N
%     Sum = Sum + 0.5 * V(j1, j2, 2) .* (Hp(k2, j2, j1) .* phi_A_(j1 - 1, j2 - 1) + Hm(k2, j2, j1) .* phi_A_(j1 - 1, -j2 -1));
% end
% Sum