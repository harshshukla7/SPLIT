%% Infinite horizon MPC using ADMM - main

clear all;
close all;
randn('state',1);
rand('state',1);

MAXITER = 10000;

%% Defining the problem parameters
n=7;
A = randn(2*n,n); b = randn(2*n,1);   c = randn(n,1); d = 2;

H = eye(n);
sigma_G = min(eig(H));

%% Solve with CVX
cvx_begin
cvx_solver 'sedumi'
    variables Z(n,1) T;
    obj = .5*Z'*H*Z;
    minimize ( obj )
    subject to 
    { A*Z + b,c'*Z + d } == lorentz(2*n);
cvx_end

%% Compare to Yalmip
yal_z = sdpvar(n, 1);
con = [];
con = set(cone(A*yal_z + b,c'*yal_z + d));
yal_obj = .5*yal_z'*H*yal_z;

result = solvesdp(con, yal_obj, sdpsettings('verbose', 0, 'solver', 'cplex'));

if result.problem, error('YALMIP could not solve problem'); end


splitProb.clearProblem;
s = splitvar(n,1);
obj = .5*s'*H*s;
norm(A*s + b,2) <= c'*s + d;
minimize(obj);
prob = splitProb.genProblem;
settings.dualTol   = 1e-5; 
settings.primalTol = 1e-5;
settings.maxItr    = 1e4;
settings.method = 'admm';
[sol_admm, stats_admm, data_to_export_admm] = admm(prob, settings);
    

    
% %     %% solve with Chambolle I
% %     lambda1 = zeros(2*n,1); old.lambda1 = lambda1; tilde.lambda1 = lambda1;
% %     lambda2 = 0; old.lambda2 = lambda2; tilde.lambda2 = lambda2;
% %     z = zeros(n,1);   old.z = z;  bar.z = z;
% %     rho = 1;
% %     
% %     
% %     %%admm
% %      v=zeros(2*n,1); t=0;
% %         for jjj = 1:2000
% %     z = (H+rho*(K'*K)) \ (K'*[lambda1;lambda2]+rho*K'*([v;t]-k));
% %   
% %     tilde.v = A*z+b-lambda1/rho;
% %     tilde.t = c'*z+d-lambda2/rho;
% % % %     u = proj_secondOrderCone([tilde.t;tilde.v]);
% % % %     t = u(1,1); v = u(2:end,1);
% %     ny = norm(tilde.v);
% % 
% %     if ny <= tilde.t
% %       v = tilde.v;
% %       t = tilde.t;
% %     elseif ny <= -tilde.t
% %       v = zeros(length(tilde.v), 1);
% %       t = zeros(length(tilde.t), 1);
% %     else
% %       v = tilde.v * (tilde.t+ny)/(2*ny);
% %       t = ny * (tilde.t+ny)/(2*ny);
% %     end
% %     
% %           
% %           
% %     lambda1 = lambda1 + rho*(v - A*z - b);
% %     lambda2 = lambda2 + rho*(t - c'*z - d);
% %     end
        
        
%         return;
    
        

%% Compute chambolle parameters
K = [A; c'];
k = [b; d];
Mu = max(svds(K)); % Bound on K operator
no.con.all = size(K,1);

    theta = 1; sigma = 1*.9/Mu; tau = (1/1)*sigma;
%     LH = chol(H+(1/tau),'lower');
    lambda1 = zeros(2*n,1); old.lambda1 = lambda1; tilde.lambda1 = lambda1;
    lambda2 = 0; old.lambda2 = lambda2; tilde.lambda2 = lambda2;
    z = zeros(n,1);   old.z = z;  bar.z = z;
    
    % reformulation of the problem with new params introduced by myself
    
    lambda1 = zeros(2*n, 1);
    lambda2 = 0;
    lam_main = [zeros(2*n, 1); 0];
    lam_prev = size(lam_main);
    lam_tilde = lam_main;
    
    % min_x max_{lambda1,lambda2} -<Az+b,lambda1> - <c'z+d,lambda2> + (1/2)z'Hz - I(|lambda_1|_2<=lambda_2)
    
    tic
    tilde_lambda = [];
    for jjj = 1:MAXITER
        
        % step 1: Solve argmin I(|lambda_1|_2<=lambda_2) + <Az+b,lambda1> + <c'z+d,lambda2> + (1/2/sigma)|(lambda_1,lambda_2)-(lambda_1,lambda_2)^k|_2^2
        lam_prev    = lam_main;
        lam_reverse = [lambda2; lambda1];
        lam_tilde = lam_reverse - sigma*(full(prob.prox.L)*z+full(prob.prox.l));
               
        % u = proj_secondOrderCone([tilde.lambda2; tilde.lambda1]);
        u = proj_secondOrderConeConj(lam_tilde);
        lambda2 = u(1,1); lambda1 = u(2:end,1);
        lam_main = [lambda1; lambda2];
                
        % step 2: Solve argmin (1/2)z'Hz - <Az+b,lambda1> - <c'z+d,lambda2> + (1/2/tau)|z-z^k|_2^2
        old.z = z;
        z = (H+eye(size(H,1))*(1/tau)) \ ((1/tau)*z+(K'*lam_main));
        
        % step 3
        bar.z = z + theta*(z-old.z);
        
        % step 4: compute residual
        res.s  = tau \ ( z-old.z);
        res.r  = K*(bar.z-old.z) - sigma \ (lam_main - lam_prev);
        rDual(jjj) = norm(res.s);    rPrimal(jjj) = norm(res.r);
        
        % check convergence 
        if (rDual(jjj) <= 1e-5 && rPrimal(jjj) <= 1e-5)
                break;
        end
    end
    toc
    
    number_of_iterations = jjj;
% %     semilogy(rPrimal,'r'); semilogy(rDual,'g');
% %     
% %     clear rDual;    clear rPrimal;

  