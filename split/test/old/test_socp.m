%% Infinite horizon MPC using ADMM - main

clear all;
close all;
randn('state',0);
rand('state',0);

% % MAXITER = 100;

%% Defining the problem parameters
n=7;
A = randn(2*n,n); b = randn(2*n,1);   c = randn(n,1); d = 2;

% % % %% Compute chambolle parameters
% % K = [A; c'];
% % k = [b; d];
% % Mu = max(svds(K)); % Bound on K operator
% % no.con.all = size(K,1);


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

%%

splitProb.clearProblem;
z = splitvar(n,1);
obj = .5*z'*H*z;
norm(A*z + b,2) <= c'*z + d;
minimize(obj);
prob = splitProb.genProblem;
settings.dualTol   = 1e-5; 
settings.primalTol = 1e-5;
settings.maxItr    = 1e4;
settings.method = 'admm';
[sol_admm, stats_admm, data_to_export_admm] = admm(prob, settings);
    

% % % %     %% solve with Chambolle I
% % % %     lambda1 = zeros(2*n,1); old.lambda1 = lambda1; tilde.lambda1 = lambda1;
% % % %     lambda2 = 0; old.lambda2 = lambda2; tilde.lambda2 = lambda2;
% % % %     z = zeros(n,1);   old.z = z;  bar.z = z;
% % % %     
% % % %     
% % % %     %%admm
% % % %     rho=1; v=zeros(2*n,1); t=0;
% % % %         for k = 1:1000
% % % %     z = (H+rho*(K'*K)) \ (rho*K'*[lambda1;lambda2]+rho*K'*([v;t]+k));
% % % %     tilde.v = A*z+b+lambda1;
% % % %     tilde.t = c'*z+d+lambda2;
% % % %     u = proj_secondOrderCone([tilde.t;tilde.v]);
% % % %     t = u(1,1); v = u(2:end,1);
% % % % % %             if norm(tilde.v) <= -tilde.t
% % % % % %             v = 0;t = 0;
% % % % % %             elseif norm(tilde.v) <= tilde.t
% % % % % %             v= tilde.v; t = tilde.t;
% % % % % %             elseif ( norm(tilde.v)>abs(tilde.t) )
% % % % % %             v = .5*(1+(tilde.t/norm(tilde.v))) * tilde.v;
% % % % % %             t = .5*(1+(tilde.t/norm(tilde.v))) * norm(tilde.v);
% % % % % %             end
% % % %         lambda1 = tilde.v-v;
% % % %         lambda2 = tilde.t-t;
% % % %         end
% % % %         
% % % %         
% % % %         return;
% % % %     
% % % %         
% % % %     theta = 1; sigma = 1*1/Mu; tau = (1/1)*sigma;
% % % %     LH = chol(H+(1/tau),'lower');
% % % %     lambda1 = zeros(2*n,1); old.lambda1 = lambda1; tilde.lambda1 = lambda1;
% % % %     lambda2 = 0; old.lambda2 = lambda2; tilde.lambda2 = lambda2;
% % % %     z = zeros(n,1);   old.z = z;  bar.z = z;
% % % %     
% % % %     tic
% % % %     for k = 1:MAXITER 
% % % %         %% step 1
% % % %         old.lambda1 = lambda1;
% % % %         old.lambda2 = lambda2;
% % % %         tilde.lambda1 = lambda1 + sigma*(A*z+b);
% % % %         tilde.lambda2 = lambda2 + sigma*(c'*z+d);
% % % % %         if norm(tilde.lambda1) <= -tilde.lambda2
% % % % %             lambda1 = tilde.lambda1; lambda2 = tilde.lambda2;
% % % % %         elseif ( norm(tilde.lambda1)>abs(-tilde.lambda2) )
% % % % %             lambda1= .5*(1-tilde.lambda2/norm(tilde.lambda1)) * tilde.lambda1;
% % % % %             lambda2 = .5*(1-tilde.lambda2/norm(tilde.lambda1)) * norm(tilde.lambda1);
% % % % %         elseif norm(tilde.lambda1) <= tilde.lambda2
% % % % %             lambda1 = 0; lambda2 = 0;
% % % % %         end
% % % %     cvx_begin
% % % %     cvx_solver 'sdpt3'
% % % %         variables L1(2*n) L2(1)
% % % %         obj = .5*quad_form([L1;L2]-[tilde.lambda1;tilde.lambda2],eye(2*n+1));
% % % %         minimize ( obj )
% % % %         subject to 
% % % %          norm(L1,2) <= -L2;
% % % %     cvx_end
% % % %         lambda1=L1; lambda2=L2;    
% % % %         
% % % %         %% step 2
% % % %         old.z = z;
% % % %         z = LH' \ (LH \ ((1/tau)*z-(K'*[lambda1;lambda2])) );
% % % %         
% % % %         %% step 3
% % % %         bar.z = z + theta*(z-old.z);
% % % %         
% % % %         res.s  = tau \ ( z-old.z);
% % % %         res.r  = K*(bar.z-old.z) - sigma \ ([lambda1;lambda2]-[old.lambda1;old.lambda2]);
% % % %         rDual(k) = norm(res.s);    rPrimal(k) = norm(res.r);
% % % %         if (rDual(k) <= 1e-3 && rPrimal(k) <= 1e-3)
% % % %                 break;
% % % %         end
% % % %     end
% % % %     toc
% % % %     
% % % %     pause;
% % % %     k
% % % %     hold;
% % % %     semilogy(rPrimal,'r'); semilogy(rDual,'g');
% % % %     
% % % %     clear rDual;    clear rPrimal;

  