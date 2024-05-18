function [yh, yl] = random_walk(N,true_vol,high_stc, low_stc)
x = zeros(N,1);

for t=2:N
    x(t) = x(t-1) + sqrt(true_vol)*randn;   % adds volatility to reward rate
end

yh = x + sqrt(high_stc)*randn(N,1);
yl = x + sqrt(low_stc)*randn(N,1);      % adds stochasticity to outcome

end
