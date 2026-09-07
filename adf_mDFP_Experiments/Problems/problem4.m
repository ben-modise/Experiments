function y = problem4(x)
n = length(x);
y = 2*x +sin(x) -1 -[x(2:n);0];