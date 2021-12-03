function M = get_M(x1, x2, a, b, N)
%{
    Description: Computes a matrix containing all H(k2, j2, j1) values. 
    
    Parameters:
      x1:       [1x1 real] Lower integration range
      x2:       [1x1 real] Upper integration range
      a:        [1x1 real] Cosine argument
      b:        [1x1 real] Cosine argument
      N:        [1x1 integer] Truncation element


    Output: 
      M:        [N x N] matrix containing all M's 
    
    References:

%}
pbma = pi / b - a;
[k, j] = deal(0:(N-1));

%% Compute MC
jpk = k'+ j; %Compute in vectorized form
Mc = 1 ./ jpk .* (exp(1i * jpk * (x2 - a) * pbma) ...
        - exp(1i * jpk * (x1 - a) * pbma));
Mc(1,1) = 1i * pbma * (x2 - x1); %Correct NA terms

%% Compute MS
jmk = k'- j; %Compute in vectorized form
Ms = 1 ./ jmk .* (exp(1i * jmk * (x2 - a) * pbma) ...
        - exp(1i * jmk * (x1 - a) * pbma));
for idx = 1:length(k)
    Ms(idx, idx) = 1i * pbma * (x2 - x1); %Correct NA terms
end

%% Compute M
M = - 1i / pi * (Mc + Ms); %return the result

end