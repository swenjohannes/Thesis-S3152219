function y = bsprice(S,K,r,c,vol,T,pc)

vsqrt = vol * sqrt(T + 1e-20);
d1 = (log(S./K) + (r-c)*T)/vsqrt +vsqrt/2;
d2 = d1 - vsqrt;

y = pc * (exp(-c*T) * S .* normcdf(d1*pc,0,1)...
  - exp(-r*T) * K .* normcdf(d2*pc,0,1));
end
