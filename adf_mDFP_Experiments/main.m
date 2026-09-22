% Runs test for adf-mDFP
% 
addpath('Problems')
n = 1000;
dimensions = [10000, 50000, 100000];
problems = {@problem1, @problem2, @problem3, @problem4...
    @problem5, @problem6, @problem7, @problem8};
%Starting points here use Case III
%[x0, x_m1] = generateStartingPoints(1,n);


%maxIter = 100;
%tolerance = 1e-6;
%J = problems{4};

%[x_star1, k1, fun_evals1, Jnorm1, iflag1] = BNN(J, x0);

output = fopen("results.txt", "wt");
fprintf(output, '\n Experiment results \n\n');
fprintf(output, '\n ---------------------------------------------------------\n');
fprintf(output, 'Case \t\t ObjFun\t\t Dim\t iflag \t NI \t NFE\t\t CPU \t\t normF\n');

for prob = 1:8
    J = problems{prob};
    for i = 1:7
        Case = intToRoman(i);
        [x0, x_m1] = generateStartingPoints(i, n);
        fprintf(output, '%s\t\t\t%9s\t%5d\t', Case, func2str(J), n);
    
        tic;
        [x_star, k, fun_evals, Jnorm, iflag] = adf_mDFP(J, x0, x_m1);
        tsec = toc;
        fprintf(output,'%2d\t\t %3d\t%4d\t%10.8f\t\t%12.8f\n',iflag,k,fun_evals,tsec,Jnorm);
    end
end


fclose('all');