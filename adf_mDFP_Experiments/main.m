% Runs test for adf-mDFP
% 
addpath('Problems\')
n = 1000;
dimensions = [10000, 50000, 100000];
%Starting points here use Case III
x0 = 0.5*ones(n,1);
x_m1 = 0.1*ones(n,1);

%maxIter = 100;
%tolerance = 1e-6;
J = @problem5;
%[x_star, k, fun_evals, Jnorm] = adf_mDFP(J, x_m1, x0);
