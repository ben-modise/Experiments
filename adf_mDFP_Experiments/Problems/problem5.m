function y = problem5(x)

n = length(x);
m = diag(5/2*ones(n,1)); %main diagonal
a = diag(ones(n-1,1),1);
b = diag(ones(n-1,1),-1);
A = m+a+b;

y = A*x - 1;