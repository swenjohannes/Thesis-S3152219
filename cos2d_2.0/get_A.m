function A = get_A(Hp, Hm, Wp, Wm, N)
%{
    Description: Computes the A matrix 
    
    Parameters:
      Hp:       [j1 x j2 x k2 complex] Positive part of H matrix
      Hm:       [j1 x j2 x k2 complex] Negative part of H matrix
      Wp:       [j1 x j2 x k2 complex] Positive part of H matrix
      Wm:       [j1 x j2 x k2 complex] Negative part of H matrix
      N:        [1x1 real] Truncation argument

    Output: 
      A:        [j1 x k2] matrix containing weighted_sum(Hp * Wp + Hm * Wm)
    
    References:

%}

Wp(:,1) = 0.5 * Wp(:,1);  % Set first j2 element of both W matrices to 
Wm(:,1) = 0.5 * Wm(:,1);  % a half to compute the weighted sum!
A = reshape(sum(Hp .* Wp + Hm .* Wm, 2), [N N]); %Take rowsums and remove 3D!
end
