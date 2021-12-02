function Sum = weighted_sum(x)
    Sum = 0.5 * x(1) + sum(x(2:end)); %Weight first term by a half
end