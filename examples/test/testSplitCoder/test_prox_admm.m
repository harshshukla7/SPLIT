clear all

n = 5;
m = 5;
N = 15;

Q = eye(n);
R = eye(m);

F = ones(n*N,n);
q = ones(n*N,1);
f = randn(n,1);
d = 1;

K = randn(3*(n*N+m*(N-1)),n*N+m*(N-1));
k = abs(randn(3*(n*N+m*(N-1)),1));

[A,B,C,D] = ssdata(drss(n,m,m));

x = splitvar(n,N);
u = splitvar(m,N-1);


%% 

clear all

n = 5;
m = 5;
N = 15;

[A,B,C,D] = ssdata(drss(n,n,m));

x = splitvar(n,N);
u = splitvar(m,N-1);
x0 = parameter(n,1);

x(:,1) == x0;


%%%
F = ones(n*N,n);
q = ones(n*N,1);
f = randn(n,1);
K = randn(3*(n*N+m*(N-1)),n*N+m*(N-1));
k = abs(randn(3*(n*N+m*(N-1)),1));
K*[x(:); u(:)] <= k;

%%-10 <= x(:) <= 10;
%%-1 <= u(:) <= 1;

x(:,2:end) == A*x(:,1:end-1) + B*u;

Q = eye(n);
R = eye(m);
objX = x .* (Q*x);
objU = u .* (R*u);
obj = sum(objX(:)) + sum(objU(:));

minimize(obj);

prob = splitProb.genProblem;

%% KKT Matrix

Q = ones(n,n);



%% box ... very sparse
-1 <= x(:) <= 1;
-1 <= u(:) <= 1;
x(:,2:end) == A*x(:,1:end-1) + B*u;

minimize(trace(x'*Q*x) + trace(u'*R*u))
clear prob
prob = splitProb.genProblem;
figure()
spy(prob.prox.L)

figure()
spy(prob.dat.Q)
title('KKT')
%% polytope ... fully populated
K*[x(:); u(:)] <= k;
x(:,2:end) == A*x(:,1:end-1) + B*u;

objX = x .* (Q*x);
objU = u .* (R*u);
obj = sum(objX(:)) + sum(objU(:));

minimize(obj);


prob = splitProb.genProblem;
figure()
spy(prob.prox.L)


cd = coder_admm(prob);
%% Lorentz cone .. very sparse
norm(F'*x(:)+f,2) <= q'*x(:) + d;
x(:,2:end) == A*x(:,1:end-1) + B*u;
minimize(trace(x'*Q*x) + trace(u'*R*u))
clear prob
prob = splitProb.genProblem;
figure()
spy(prob.prox.L)

d = coder_admm(prob);

%% l1-ball .. very sparse exactly like above
norm(F'*x(:)+f,1) <= d;
x(:,2:end) == A*x(:,1:end-1) + B*u;

minimize(trace(x'*Q*x) + trace(u'*R*u))
clear prob
prob = splitProb.genProblem;
figure()
spy(prob.prox.L)

d = coder_admm(prob);
%% l2-ball the same as above ! sparse
norm(F'*x(:)+f,2) <= d;
x(:,2:end) == A*x(:,1:end-1) + B*u;

minimize(trace(x'*Q*x) + trace(u'*R*u))
clear prob
prob = splitProb.genProblem;
figure()
spy(prob.prox.L)

sd = coder_admm(prob);

%% linfty-ball again the same as above sparse ! 
norm(F'*x(:)+f,Inf) <= d;
x(:,2:end) == A*x(:,1:end-1) + B*u;

minimize(trace(x'*Q*x) + trace(u'*R*u))
clear prob
prob = splitProb.genProblem;
figure()
spy(prob.prox.L)

d = coder_admm(prob);
%% equality constraints
x(:,2:end) == A*x(:,1:end-1) + B*u;

minimize(trace(x'*Q*x) + trace(u'*R*u))
clear prob
prob = splitProb.genProblem;
figure()
spy(prob.prox.L)

% % %%
% % fcnList_target = { @fama @admm @admm_prec @fama_prec };
% % 
% %     settings.maxItr          = 5e4;
% %     settings.dualTol         = 1e-3;
% %     settings.primalTol       = 1e-3;
% %     settings.fun_eval        = 1;
% % 
% %     if ( kkk ==  2 )
% %         [rho, alpha, P2] = optimal_selections_linear_ADMM(K,eye(nxu),k);
% %         settings.rho = rho;
% %         settings.relaxation = 1.9;
% %     elseif ( kkk ==  3 )
% %         
% %         [rho, alpha, P2] = optimal_selections_linear_ADMM(K,eye(nxu),k);
% %         settings.rho = rho;
% %         settings.relaxation = 1.9;
% %     elseif ( kkk == 4 )
% %         settings.rho             = 1;
% %         settings.dualTol         = 1e-3;
% %         settings.primalTol       = 1e-3;
% %         settings.scale_rPrimal   = eye(size(L',1));
% %     end
% %     
% % [sol, stats, stats_plot, ~] = fcnList_target{kkk}(prob, settings);
% % sol.x = (L') \ sol.x;
% % SOL{kkk,ttt} = sol;
% % ITER_COUNT(kkk,ttt) = stats.numiter;
% % STATS_PLOT{kkk,ttt,:} = stats_plot;
% % fom.x_struct{kkk,ttt} = sol.x(1:SimpleModel.nx);
% % fom.u_struct{kkk,ttt} = sol.x(SimpleModel.nx+1:end);
% % fom.x{kkk} = [[fom.x{kkk}] [fom.x_struct{kkk,ttt}]];
% % fom.u{kkk} = [[fom.u{kkk}] [fom.u_struct{kkk,ttt}]];
% % 
