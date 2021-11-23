function mj = mj(a, b, x1, x2, j)
%{
 Inputs : a             - Cosine argument (lower truncation bound) 
        : b             - Cosine argument (upper truncation bound)
        : x1             
        : x2             
        : j             

%}
pbma = pi / (b - a);
u = (x2 - a) * pbma;
v = (x1 - a) * pbma;

mj = (exp(1i * j * u) - exp(1i * j * v)) ./ j;
mj(j == 0) = (x2 - x1) * pbma * 1i;