function [x,k,numf,gnorm, iflag] = BD(f,x0,P,P_options,params)
%   Hybrd Conjugate Gradient Projection method for convex constrained monotone
%   systems. Based on 'A modified Hestenes–Stiefel conjugate gradient method with
%   sufficient descent condition and conjugacy condition' by Dong et al(2015).
%   
%
%   This program implements the algorithm for the above paper for 
%   solving convex constrained equations that are monotone.

%   Input: x0  = initial iterate -- must be within the constraints
%          f   = objective function f: R^n-->R^n
%          P   = the projection function onto the constraint region C
%                optional, default: P(x) = x
%    P_options = any optional arguments to be used in the projection
%                function
%          tol = termination criterion norm(F_k) <= tol
%                optional, default = 1.e-5
%        maxit = maximum iterations (optional) default = 50000
%
%   Output: x     = solution
%           k  = iteration history
%           numf  = number of function evaluations
%           norm(Fk) = norm of objective function value
%           iflag = 0 algorithm successful
%                 = 1 steplength unsuccessful
%                 = 2 norm(Fk) > tol && iter > maxit
    
    
if nargin < 3
    P = @constantfunction;
elseif isempty(P)
    P = @constantfunction;
end
keyboard;
if nargin >= 4 && ~isempty(P_options)
    P = @(x)P(x,P_options);
end

% input validation

if norm(P(x0)-x0)~=0
    error('x0 is not in C');
end


% initialize variables
if nargin < 5
    maxit = 2000;
    tol = 1.e-5;
elseif isempty(params)
    maxit = 2000;
    tol = 1.e-5;
else
    maxit = params(1);
    tol = params(2);
end


%TODO: perhaps the following should also be inputs with default values

rho = 0.5;
sigma = 1e-4;
eta = 1.0;
ga = 1.85; %gamma variable
%b_min = 1;
%b_max = 10;
k = 0;
bk = 1; %beta for the line search (whereas b_bd is for computing direction)
t = 0.1;

xk = x0; 
Fk = feval(f,xk);
numf = 1;
iflag = 0;
% variablename_1 represents the value at k-1
Fk_1 = Fk;
dk_1 = zeros(size(Fk));
%alphak_1 = 0;

%b_bd = 0;
%A = 1; %used to store Fk'*y_k-1 value
B = 1; %used to store Fk'*d_k-1 value

count1 = -1;
count2 = 0;
count3 = 0;

while(norm(Fk) > tol && k <= maxit) %step1
    % step2: calculation of direction
    if isequal(k,0) || A < 0
        %b_bd=0;
        dk = -Fk;
        count1 = count1 + 1;
    elseif B > 0  
        C = dk_1'*yk_1;   
        %calculation of the hybrid beta parameter
        etak_1 = -1/(norm(dk_1)*min([eta, norm(Fk_1)]));
        b_hs = (Ab)/(C);
        b_dhs = (-(Fk_1'*dk_1)/(C))*b_hs - t*(norm(ybk_1)^2*B)/(C^2);
        b_D = max([etak_1,b_dhs]);
        lambda_k = 1 + (B*(Fk'*yk_1))/(C*norm(Fk)^2);
        dk = -lambda_k*Fk + b_D*dk_1;
        count2 = count2 + 1;
    elseif B <= 0
        
        %proved using Fk'*yk_1
        bb_hs = (Fk'*yk_1)/(dk_1'*ybk_1);
        dk = -Fk + bb_hs*dk_1;
        count3= count3 + 1;
        
    end

    % we cannot use this if we do not have that dk is bounded above
    % by a multiple of norm(Fk)
    k = k+1;
    if norm(dk)<tol
        x = xk; 
        gnorm = norm(Fk);
        return;
    end
    
    % step3: finding the trial point zk - line search
    mk=0;
    alpha_k = bk*rho^mk;
    zk = xk + alpha_k*dk;  
    Fz = feval(f,zk);
    numf = numf+1;
    %varsigma=1;
    
    while(-Fz'*dk < sigma*alpha_k*norm(dk)^2) %*varsigma
        mk = mk+1;
        if(mk > 1000)
            disp('trial point zk iterations exceeded')
            iflag = 1;
            return;
        end
        alpha_k = bk*rho^mk;
        zk = xk + alpha_k*dk;  
        Fz = feval(f,zk);
        %varsigma = norm(Fz)/(1+norm(Fz));
        numf = numf+1;
    end
    %k = k+1;
    %check if Fz==0 and stops if true
    if norm(Fz)<1.e-6  %has to be tighter than tol since Fk would not have been updated
        x = xk; 
        gnorm = norm(Fk);
        return
    end
    
    % step4: computing the next iterate 
    Xi = (Fz'*(xk - zk))/norm(Fz)^2;
    xk = P(xk - ga*Xi*Fz);
    
    % Set bk to be in [bmin, bmax]. For now we let bk=1 for all k
    
    alphak_1 = alpha_k;
    Fk_1 = Fk;
    dk_1 = dk;
    Fk = feval(f, xk);
    numf = numf+1;
    yk_1 = Fk - Fk_1;
    deltak_1 = 1 + norm(Fk_1)^-1*max([0,((-alphak_1*dk_1'*yk_1)/(alphak_1^2*norm(dk_1)^2))]);
    ybk_1 = yk_1 + deltak_1*alphak_1*norm(Fk_1)*dk_1;
    A = Fk'*yk_1;
    Ab = Fk'*ybk_1;
    B = Fk'*dk_1;
    
end
%disp([count1,count2,count3])
x = xk; 
gnorm = norm(Fk);
%k = k-1;
if (norm(Fk) > tol && k >= maxit)
    iflag =2;
end
