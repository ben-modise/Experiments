function y = problem6(x)

n = length(x);
%y = zeros(n, 1); % Initialize output array
y = exp(x) -1 + [0;x(2:n)];