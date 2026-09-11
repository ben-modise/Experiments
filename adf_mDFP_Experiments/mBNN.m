function [x,k,numf,gnorm, iflag] = mBNN(f,x0)
% Modified BNN


rho = 0.5;
sigma = 1e-4;
mu = 1.4;
ga = 1.85; %gamma variable
%b_min = 1;
%b_max = 10;
k = 0;
bk = 1; %beta for the line search (whereas b_nn is for computing direction)

xk = x0; 
Fk = feval(f,xk);
numf = 1;
iflag = 0;
% variablename_1 represents the value at k-1
Fk_1 = Fk;
dk_1 = zeros(size(Fk));

%b_nn = 0;

while(norm(Fk) > tol && k <= maxit) %step1
    % step2: calculation of direction
    if k==0
        b_nn=0;
    else
        Fk = feval(f,xk);
        numf = numf + 1;
        %calculation of the hybrid beta parameter
        A = norm(Fk)^2;
        temp1 = (norm(Fk)/norm(Fk_1))*(Fk'*Fk_1);
        B = max([0,temp1]);
        temp2 = mu*norm(dk_1)*norm(Fk);
        temp3 = dk_1'*(Fk-Fk_1);
        C = max([temp2,temp3]);
        b_nn = (A-B)/C;
    end
    dk = -Fk + b_nn*dk_1;

    if norm(dk)<tol
        x = xk; 
        gnorm = norm(Fk);
        return;
    end

    % step3: finding the trial point zk - line search
    mk=0;
    alp = bk*rho^mk;
    zk = xk + alp*dk;  
    Fz = feval(f,zk);
    numf = numf+1;
    while(-Fz'*dk < sigma*alp*norm(dk)^2)
        mk = mk+1;
        if(mk > 1000)
            disp('trial point zk iterations exceeded')
            iflag = 1;
            return;
        end
        alp = bk*rho^mk;
        zk = xk + alp*dk;  
        Fz = feval(f,zk);
        numf = numf+1;
    end
    %the algorithm given checks if Fz==0 and stops if true
    if norm(Fz)<1.e-6  %has to be tighter than tol since Fk would not have been updated
        x = xk; 
        gnorm = norm(Fk);
        return
    end
    k = k+1;
    % step4: computing the next iterate 
    Xi = (Fz'*(xk - zk))/norm(Fz)^2;
    xk = P(xk - ga*Xi*Fz);

    % Set bk to be in [bmin, bmax]. For now we let bk=1 for all k

    Fk_1 = Fk;
    dk_1 = dk;
end
x = xk; 
gnorm = norm(Fk);
%k = k-1;
if (norm(Fk) > tol && k >= maxit)
    iflag =2;
end