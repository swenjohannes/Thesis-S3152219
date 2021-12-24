function x = bsimpvol(S,K,r,c,T,pc,price)

x = zeros(size(K));
for ii = 1:length(K)
    myfun = @(v) bsprice(S,K(ii),r,c,v,T,pc)-price(ii);
    x(ii) = fzero(myfun, 0.1);
end
