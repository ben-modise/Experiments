function y = problem2(x)
n = length(x);
y = zeros(n, 1); % Initialize output array
for i = 2:n-1
    y(i) = x(i) - exp(cos((x(i-1) + x(i)+x(i+1))/(i)));
end
y(1) = x(1) - exp(cos((x(1)+x(2))/(2)));
y(n) = x(n) - exp(cos((x(n)+x(n-1))/(n)));