function [x_s, k, fevals, normJ] = adf_mDFP(J, x_prev, x0, ~)
%=========================================================================
% adf_mDFP
%
% Accelerated derivative-free memoryless DFP method
%
% Implements proposed algorithm in
%
%"An accelerated derivative-free memoryless Davidon-Fletcher-Powell 
% method and its iteration-complexity analysis."
%  
% Lui, Li, Shao, and Wu (2025)
%==========================================================================

%% ------------- Initialize parameters---------------------------------
% TODO: Set using structure using function defaultOptions()
constants = {0.01, 1, 0.6, 1.7, 0.01, 0.001, 0.8, 0.001};
[sigma, alpha, tao, rho, varrho, mu1, mu2, r] = constants{:};

tol1 = 1e-7; % for checking against norm(dk)
tol2 = 1e-6; % for checking against norm(Jk)
maxit = 10000; %; % maximum number of iterations
maxLineSearchTrials = 10000;



%%-----------------------------------------------------------------------\
%% --------------------subfunctions--------------------------------------
% Necessary functions for different steps in the main iteration
% lineSearch, searchDirection


    function proj_nu = P(nu)
        %implements box projection for scalar nu
        %between mu1 and mu2
        proj_nu = max(mu1, min(nu, mu2)); % Box projection for nu
    end

    function [zk, J_zk] = lineSearch()
        % Perform line search to find optimal step size
        % assumes dk, norm_dk, omega and J exist in variable space
        alpha_k = alpha; % Initial step size
        norm_dk2 = norm_dk^2;
        zk = omega + alpha_k*dk;
        J_zk = J(zk); fevals = fevals + 1;
        lineSearchtrials = 0;
        while -J_zk'*dk < sigma*alpha_k*P(norm(J_zk))*norm_dk2
            alpha_k = alpha_k * tao; % Update step size
            zk = omega + alpha_k*dk;
            J_zk = J(zk); fevals = fevals + 1;
            lineSearchtrials = lineSearchtrials +1;
            if lineSearchtrials >= maxLineSearchTrials
                disp('line search trials exceeded')
                break
            end
        end
    end

    function dk = computeDirection()        
        if k>=0
            dk = -Jomega;
            return
        end
        sk_1 = omega - omega_prev;
        lambda_option = -((Jomega - Jomega_prev)'*sk_1)/norm(sk_1)^2;
        lambda = r + max(0, lambda_option);
        yHat = Jomega - Jomega_prev + lambda*sk_1;

        if norm(yHat) < tol2
            dk = -Jomega;
            return
        end
        sk_yHat = sk_1'*yHat;
        psiOptionB = norm(sk_1)^2/(sk_yHat);
        psiOption = norm(sk_1)^4*norm(yHat)^2;
        denom = 2*(sk_yHat)^3 - norm(sk_1)^2*norm(yHat)^2*(sk_yHat);
        psiOption = psiOption/denom;
        psi = max(psiOption, psiOptionB);
        term1 = -Jomega;
        term2 = ((Jomega'*yHat)/norm(yHat)^2)*yHat;
        term3 = - (1/psi)*(Jomega'*sk_1)/(sk_yHat)*sk_1;
        dk = term1 + term2 + term3;
        % if k==1
             % keyboard;
        % end
    end
%% -----------------------------------------------------------------------

%tic
%% ==============Main Iteration===========================================
xk = x0;
fevals = 0;
for k = 0:maxit
    % ------------------------------------------------
    % Step 1: Inertial extrapolation       (Eq 14)
    % ------------------------------------------------
    omega = xk + varrho*(xk - x_prev);
    Jomega = J(omega); fevals = fevals+1;

    if norm(Jomega) < tol2
        x_s = omega;
        normJ = norm(Jomega);
        %keyboard;
        break
    end

    % ------------------------------------------------
    % Step 2: Search direction            (Eq 15)
    % ------------------------------------------------

    % s = omega - omega_prev;
    dk = computeDirection();
    %TODO: Debug, check if direction satisfies sufficient descent
    norm_dk = norm(dk);
    if norm_dk < tol1
        x_s = omega;
        normJ = norm(Jomega);
        %keyboard;
        break
    end
    % -------------------------------------------------
    % Step 3: Line search                  (Eq 17) 
    % -------------------------------------------------
    [zk, J_zk] = lineSearch();

    if norm(J_zk) < tol2
        normJ = norm(J_zk);
        x_s = zk; 
        %keyboard; 
        break
    end

    % ------------------------------------------------------------------
    % Step 4: Hyperplane projection              (Eq 18)
    % ------------------------------------------------------------------
    %x_new = zk;
    zeta = J_zk'*(omega - zk)/norm(J_zk)^2;
    x_new = omega - rho*zeta*J_zk;

    J_xk = J(x_new);
    if norm(J_xk) < tol2
        normJ = norm(J_xk);
        x_s = x_new;
        break
    end


    % -------------------------------------------------------------------
    % Prepare next iteration
    % -------------------------------------------------------------------
    omega_prev = omega;
    Jomega_prev = Jomega;
    x_prev = xk;
    xk = x_new;

end
%% ======================================================================

if k==maxit
    normJ = norm(Jomega);
    x_s = omega;
    disp('maximum iterations reached')
end

end