% Runs test for adf-mDFP
% 
addpath('Problems\')
n = 1000;
dimensions = [10000, 50000, 100000];
problems = {@problem1, @problem2, @problem3, @problem4...
    @problem5, @problem6, @problem7, @problem8};
%Starting points here use Case III
[x0, x_m1] = generateStartingPoints(3,n);


%maxIter = 100;
%tolerance = 1e-6;
J = problems{7};
[x_star, k, fun_evals, Jnorm] = adf_mDFP(J, x_m1, x0);
