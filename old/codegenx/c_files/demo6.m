clear

clear all;
close all;
randn('state',1);
rand('state',0);

% Must always call clearProblem before starting! 
% (or use *clear all*, clear is not enough)
splitProb.clearProblem;

N = 6;
n = 5;
p = 2;
m = 3;
[A,B,C,D] = ssdata(drss(n,p,m));

% A = rand(n,n);
% B = rand(n,m);

x = splitvar(n, N);
u = splitvar(m, N-1);

% that would be the case, once we implement inverse in C
% via soft proposed by Giorgos/Jean
% Q = randn(n); Q = Q*Q';
% R = randn(m); R = R*R';

% now we just consider diagonal matrices in the objective
% result of the assumption: simpler C implementation
% i.e. instead of inverting a matrix one has to compute 
% only a simple addition between two vectors

% generate a random vector for matrices Q and R
Qv = rand(n,1)*3.25 - 0.02*rand(n,1)*6.5;
Rv = rand(m,1)*1.25 - 0.22*rand(m,1)/0.95;

% diagonal weighting matrices for states and inputs
Q = diag(Qv);
R = diag(Rv);

obj = 0;
for i = 1:N-1
  x(:,i+1) == A*x(:,i) + B*u(:,i);
  obj = obj + x(:,i)'*Q*x(:,i) + u(:,i)'*R*u(:,i);
end
obj = obj + x(:,end)'*x(:,end);
obj = obj + 5.5*norm(x(:,2), 2);

-5  <= x <= 5;
-.5 <= u <= .5;

% % SOC constraint
norm(0.5*x(:,2) + 2, 2) + x(1,2)/4 <= 10*x(2,2);

% ellipsoidal constraint
x(:,end)'*x(:,end) <= 1e-4;

% Inf-normball constraint
norm(u(:,1), Inf) <= 0.4;

x0 = 1.5*ones(n,1);
x(:,1) == x0;

minimize(obj);

%% Generate the problem data
prob = splitProb.genProblem;

%% Solve the problem

settings.dualTol     = 1e-4; 
settings.primalTol   = 1e-4;
settings.maxItr      = 1e5;
% settings.tau         = 1;
% settings.sigma       = 1;
% settings.theta       = 1;
% settings.gamma       = min(eig(dat.Q));

% possible values: admm, ama, fadmm, fama, chambolle1 
settings.method = 'fama';

% solve admm
if strcmp(settings.method, 'admm')
    
    % solve the problem in matlab and prepare data for export
    [sol_admm, stats_admm, stats_admm_plot, data_to_export_admm] = admm(prob, settings);
    
    % export data to a header file
    export_matrices_ADMM(data_to_export_admm);
    
    %TODO: code generation
      
    if sol_admm.problem
       error('Could not solve the problem with ADMM - infeasible?')
    end
end

% solve fadmm
if strcmp(settings.method, 'fadmm')
    
    % solve the problem in matlab and prepare data for export
    [sol_fadmm, stats_fadmm, stats_fadmm_plot, data_to_export_fadmm] = fadmm(prob, settings);
    
    % export data to a header file
    export_matrices_FADMM(data_to_export_fadmm);
    
    %TODO: code generation
      
    if sol_fadmm.problem
       error('Could not solve the problem with FADMM - infeasible?')
    end 
end

% solve ama
if strcmp(settings.method, 'ama')
    
    % solve the problem in matlab and prepare data for export
    [sol_ama, stats_ama, stats_ama_plot, data_to_export_ama] = ama(prob, settings);
    
    % export data to a header file
    export_matrices_AMA(data_to_export_ama);
    
    %TODO: code generation
        
    if sol_ama.problem
       error('Could not solve the problem with AMA - infeasible?')
    end 
end

% solve fama
if strcmp(settings.method, 'fama')
    
    % solve the problem in matlab and prepare data for export
    [sol_fama, stats_fama, stats_fama_plot, data_to_export_fama] = fama(prob, settings);
    
    % export data to a header file
    export_matrices_FAMA(data_to_export_fama);
  
    %TODO: code generation
    
    if sol_fama.problem
       error('Could not solve the problem with FAMA - infeasible?')
    end 
end


% solve chambolle-pock preconditioned I
if strcmp(settings.method, 'CPIprec_G')
    
    % solve the problem in matlab and prepare data for export
    [sol_CPIprec, stats_CPIprec, data_to_export_CPIprec] = CPIprec_G(prob, settings);
    
    % export data to a header file
    export_matrices_CPIprec(data_to_export_CPIprec);
        
    if sol_CPIprec.problem
       error('Could not solve the problem with CPIprec - infeasible?')
    end 
end

% solve chambolle-pock II
if strcmp(settings.method, 'CPII_G')
    
    % solve the problem in matlab and prepare data for export
    [sol_CPII, stats_CPII, stats_CPII_plot, data_to_export_CPII] = CPII_G(prob, settings);
    
    % export data to a header file
    export_matrices_CPII(data_to_export_CPII);
        
    if sol_CPII.problem
       error('Could not solve the problem with CPII - infeasible?')
    end 
end

% solve PDHG adaptive step
if strcmp(settings.method, 'PDHG_adapt')
    
    % solve the problem in matlab and prepare data for export
    [sol_PDHG_adapt, stats_PDHG_adapt, stats_PDHG_adapt_plot, data_to_export_PDHG_adapt] = PDHG_adapt(prob, settings);
    
    % export data to a header file
    export_matrices_PDHG_adapt(data_to_export_PDHG_adapt);
        
    if sol_PDHG_adapt.problem
       error('Could not solve the problem with PDHG_adapt - infeasible?')
    end 
end

% solve PDHG backtracking
if strcmp(settings.method, 'PDHG_back')
    
    % solve the problem in matlab and prepare data for export
    [sol_PDHG_back, stats_PDHG_back, stats_PDHG_back_plot, data_to_export_PDHG_back] = PDHG_backtracking(prob, settings);
    
    % export data to a header file
    export_matrices_PDHG_back(data_to_export_PDHG_back);
        
    if sol_PDHG_back.problem
       error('Could not solve the problem with PDHG_back - infeasible?')
    end 
end

clf
% ADMM
if strcmp(settings.method, 'admm')
    semilogy([stats_admm_plot.rDual], 'b');
    hold on; grid on
    semilogy([stats_admm_plot.rPrimal], 'b-.');
    hold on; grid on
    legend('Dual residual ADMM', 'Primal residual ADMM');
end

% FADMM
if strcmp(settings.method, 'fadmm')
    semilogy([stats_fadmm_plot.rDual], 'b');
    hold on; grid on
    semilogy([stats_fadmm_plot.rPrimal], 'b-.');
    hold on; grid on
    legend('Dual residual FADMM', 'Primal residual FADMM');
end

% AMA
if strcmp(settings.method, 'ama')
    semilogy([stats_ama_plot.rDual], 'b');
    hold on; grid on
    % semilogy([stats_ama_plot.rPrimal], 'b-.');
    hold on; grid on
    % legend('Dual residual AMA', 'Primal residual AMA');
    legend('Dual residual AMA');
end

% FAMA
if strcmp(settings.method, 'fama')
    semilogy([stats_fama_plot.rDual], 'b');
    hold on; grid on
    % semilogy([stats_fama_plot.rPrimal], 'b-.');
    % legend('Dual residual FAMA', 'Primal residual FAMA');
    legend('Dual residual FAMA');
end

% CHAMBOLLE-POCKI preconditioned
if strcmp(settings.method, 'CPIprec_G')
    semilogy([stats_CPIprec_plot.rDual], 'b');
    hold on; grid on
    semilogy([stats_CPIprec_plot.rPrimal], 'b-.');
    legend('Dual residual Chambolle-Pock I precond.', 'Primal residual Chambolle-Pock I precond.');
end

% CHAMBOLLE-POCKII
if strcmp(settings.method, 'CPII_G')
    semilogy([stats_CPII_plot.rDual], 'b');
    hold on; grid on
    semilogy([stats_CPII_plot.rPrimal], 'b-.');
    legend('Dual residual Chambolle-Pock II', 'Primal residual Chambolle-Pock II');
end

% PDHG adaptive
if strcmp(settings.method, 'PDHG_adapt')
    semilogy([stats_PDHG_adapt_plot.rDual], 'b');
    hold on; grid on
    semilogy([stats_PDHG_adapt_plot.rPrimal], 'b-.');
    legend('Dual residual PDHG_adapt', 'Primal residual PDHG_adapt');
end

% PDHG backtracking
if strcmp(settings.method, 'PDHG_back')
    semilogy([stats_PDHG_back_plot.rDual], 'b');
    hold on; grid on
    semilogy([stats_PDHG_back_plot.rPrimal], 'b-.');
    legend('Dual residual PDHG_back', 'Primal residual PDHG_back');
end

%% Compare to Yalmip

yal_x = sdpvar(n, N);
yal_u = sdpvar(m, N-1);
yal_P = sdpvar(n, n, 'symmetric');

con = [];
for i = 1:N-1
  con = con + set(yal_x(:,i+1) == A*yal_x(:,i) + B*yal_u(:,i));
end

x0 = 1.5*ones(n,1);
x(:,1) == x0;

minimize(obj);

con = con + set( -5 <= yal_x <= 5) + set(-.5 <= yal_u <= .5);
con = con + set(norm(0.5*yal_x(:,2) + 2, 2) + yal_x(1,2)/4 <= 10*yal_x(2,2));
con = con + set(norm(yal_u(:,1), Inf) <= 0.4);
con = con + set(yal_x(:,end)'*yal_x(:,end) <= 1e-4);
con = con + set(yal_x(:,1) == x0);

yal_obj = 0;
for i = 1:N-1
 yal_obj = yal_obj + yal_x(:,i)'*Q*yal_x(:,i) + yal_u(:,i)'*R*yal_u(:,i);
end
yal_obj = yal_obj + yal_x(:,end)'*yal_x(:,end);
yal_obj = yal_obj + 5.5*norm(yal_x(:,2), 2);

result = solvesdp(con, yal_obj, sdpsettings('verbose', 0, 'solver', 'cplex'));

if result.problem, error('YALMIP could not solve problem'); end

%% Recover the solution

% ADMM
if strcmp(settings.method, 'admm') 
    
splitProb.setSolution(prob, sol_admm);

% Solution is now stored in the variables
cprintf('_blue', '                OPTIMAL SOLUTION from ADMM               \n')
fprintf('x = \n'); disp(full(x.val))
fprintf('u = \n'); disp(full(u.val))

cprintf('_blue', '                ERROR TO YALMIP for ADMM              \n')
fprintf('%15s = %.2e\n', '||x - yal_x||', norm(x.val-double(yal_x)));
fprintf('%15s = %.2e\n', '||u - yal_u||', norm(u.val-double(yal_u)));

end

% FADMM
if strcmp(settings.method, 'fadmm') 
   
splitProb.setSolution(prob, sol_fadmm);

% Solution is now stored in the variables
cprintf('_blue', '                OPTIMAL SOLUTION from FADMM               \n')
fprintf('x = \n'); disp(full(x.val))
fprintf('u = \n'); disp(full(u.val))

cprintf('_blue', '                ERROR TO YALMIP for FADMM              \n')
fprintf('%15s = %.2e\n', '||x - yal_x||', norm(x.val-double(yal_x)));
fprintf('%15s = %.2e\n', '||u - yal_u||', norm(u.val-double(yal_u)));

end

% AMA
if strcmp(settings.method, 'ama') 

splitProb.setSolution(prob, sol_ama);
% Solution is now stored in the variables
cprintf('_blue', '                OPTIMAL SOLUTION from AMA               \n')
fprintf('x = \n'); disp(full(x.val))
fprintf('u = \n'); disp(full(u.val))

cprintf('_blue', '                ERROR TO YALMIP for AMA               \n')
fprintf('%15s = %.2e\n', '||x - yal_x||', norm(x.val-double(yal_x)));
fprintf('%15s = %.2e\n', '||u - yal_u||', norm(u.val-double(yal_u)));

end

% FAMA
if strcmp(settings.method, 'fama') 

splitProb.setSolution(prob, sol_fama);
% Solution is now stored in the variables
cprintf('_blue', '                OPTIMAL SOLUTION from FAMA            \n')
fprintf('x = \n'); disp(full(x.val))
fprintf('u = \n'); disp(full(u.val))
 
cprintf('_blue', '                ERROR TO YALMIP for FAMA              \n')
fprintf('%15s = %.2e\n', '||x - yal_x||', norm(x.val-double(yal_x)));
fprintf('%15s = %.2e\n', '||u - yal_u||', norm(u.val-double(yal_u)));

end


% chambolle-pockI preconditioned
if strcmp(settings.method, 'CPIprec_G') 

splitProb.setSolution(prob, sol_CPIprec);
% Solution is now stored in the variables
cprintf('_blue', '                OPTIMAL SOLUTION from CHAMBOLLE-POCKI preconditioned      \n')
fprintf('x = \n'); disp(full(x.val))
fprintf('u = \n'); disp(full(u.val))
 
cprintf('_blue', '                ERROR TO YALMIP for CHAMBOLLE-POCKI preconditioned        \n')
fprintf('%15s = %.2e\n', '||x - yal_x||', norm(x.val-double(yal_x)));
fprintf('%15s = %.2e\n', '||u - yal_u||', norm(u.val-double(yal_u)));

end

% chambolle-pockII
if strcmp(settings.method, 'CPII_G') 

splitProb.setSolution(prob, sol_CPII);
% Solution is now stored in the variables
cprintf('_blue', '                OPTIMAL SOLUTION from CHAMBOLLE-POCKII      \n')
fprintf('x = \n'); disp(full(x.val))
fprintf('u = \n'); disp(full(u.val))
 
cprintf('_blue', '                ERROR TO YALMIP for CHAMBOLLE-POCKII        \n')
fprintf('%15s = %.2e\n', '||x - yal_x||', norm(x.val-double(yal_x)));
fprintf('%15s = %.2e\n', '||u - yal_u||', norm(u.val-double(yal_u)));

end

% PDHG adaptive
if strcmp(settings.method, 'PDHG_adapt') 

splitProb.setSolution(prob, sol_PDHG_adapt);
% Solution is now stored in the variables
cprintf('_blue', '                OPTIMAL SOLUTION from PDHG adaptive      \n')
fprintf('x = \n'); disp(full(x.val))
fprintf('u = \n'); disp(full(u.val))
 
cprintf('_blue', '                ERROR TO YALMIP for PDHG adaptive        \n')
fprintf('%15s = %.2e\n', '||x - yal_x||', norm(x.val-double(yal_x)));
fprintf('%15s = %.2e\n', '||u - yal_u||', norm(u.val-double(yal_u)));

end

% PDHG bakctracking
if strcmp(settings.method, 'PDHG_back') 

splitProb.setSolution(prob, sol_PDHG_back);
% Solution is now stored in the variables
cprintf('_blue', '                OPTIMAL SOLUTION from PDHG backtracking      \n')
fprintf('x = \n'); disp(full(x.val))
fprintf('u = \n'); disp(full(u.val))
 
cprintf('_blue', '                ERROR TO YALMIP for PDHG backtracking        \n')
fprintf('%15s = %.2e\n', '||x - yal_x||', norm(x.val-double(yal_x)));
fprintf('%15s = %.2e\n', '||u - yal_u||', norm(u.val-double(yal_u)));

end