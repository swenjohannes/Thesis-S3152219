function y = f(u1, u2)
for j = 1:length(u1) 
    for k = 1:length(u2)
        y(j, k) = u1(j) * u2(k);
    end
end
y = sum(y, 'all')
end