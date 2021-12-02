function M = get_M2(x1, x2, a, b, N)
bma = b - a;
pbma = pi / bma;
[k, j] = deal(0:(N-1));

jpk = k'+ j;
Mc = 1 ./ jpk .* (exp(1i * jpk * (x2 - a) * pbma) ...
        - exp(1i * jpk * (x1 - a) * pbma));
Mc(1,1) = 1i * pbma * (x2 - x1); %Correct NA terms

jmk = k'- j;
Ms = 1 ./ jmk .* (exp(1i * jmk * (x2 - a) * pbma) ...
        - exp(1i * jmk * (x1 - a) * pbma));

for idx = 1:length(k)
    Ms(idx, idx) = 1i * pbma * (x2 - x1); %Correct NA terms
end

M = - 1i / pi * (Mc + Ms);

end