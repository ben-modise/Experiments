function [x0, x_m1] = generateStartingPoints(n, dim)
%accepts cases n ranging from 1 to 7
switch n
    case 1
        x0 = 0.2*ones(dim, 1);
        x_m1 = 0.1*ones(dim, 1);
    case 2
        x0 = 0.1*ones(dim, 1);
        x_m1 = 0.2*ones(dim, 1);
    case 3
        x0 = 0.1*ones(dim, 1);
        x_m1 = 0.5*ones(dim, 1);
    case 4
        x0 = 0.1*ones(dim, 1);
        x_m1 = 1.2*ones(dim, 1);
    case 5
        x0 = 0.1*ones(dim, 1);
        x_m1 = 1.5*ones(dim, 1);
    case 6
        x0 = 0.1*ones(dim, 1);
        x_m1 = 2*ones(dim, 1);
    case 7
        x0 = 0.1*ones(dim, 1);
        x_m1 = rand(dim, 1);
end