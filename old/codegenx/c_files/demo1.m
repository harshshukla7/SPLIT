clear

clear all;
close all;

randn('state', 2);
rand('state', 2);

% Must always call clearProblem before starting! 
% (or use *clear all*, clear is not enough)
splitProb.clearProblem;

N = 7;
n = 5;
p = 2;
m = 4;

A = rand(n, n);
B = rand(n, m);

x = splitvar(n, N);
u = splitvar(m, N-1);
Q = randn(n); Q = Q*Q';
R = randn(m); R = R*R';

obj = 0;
for i = 1:N-1
  x(:,i+1) == A*x(:,i) + B*u(:,i);
  obj = obj + x(:,i)'*Q*x(:,i) + u(:,i)'*R*u(:,i);
end
obj = obj + x(:,end)'*x(:,end);
obj = obj + 10.25*norm(x(:,3), 2);

-5  <= x <= 5;
-10 <= u <= 10;

% SOC constraint
norm(0.5*x(:,2) + 2, 2) + x(1,2)/4 <= 10*x(2,2);

% Ellipsoidal constraint
x(:,end)'*x(:,end) <= 0.002;
% 
% 1-norm ball
norm(u(:,3), 1) <= 1;

x0 = 1.5*ones(n,1);
x(:,1) == x0;

minimize(obj);

%% Generate the problem data
prob = splitProb.genProblem;

%% Solve the problem

settings.dualTol    = 1e-5; 
settings.primalTol  = 1e-5;
settings.maxItr     = 1e4;
settings.relaxation = 1.5;

% possible values: admm, ama, fadmm, fama 
settings.method = 'admm';

% solve admm
if strcmp(settings.method, 'admm')
    
    % solve the problem in matlab and prepare data for export
    [sol_admm, stats_admm, stats_admm_plot, data_to_export_admm] = admm(prob, settings);
    
    % export data to a header file
    export_matrices_ADMM(data_to_export_admm);
    
    % TODO: code generation
      
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
    
    % TODO: code generation
    
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
    
    % TODO: code generation
    
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
  
    % TODO: code generation
    
    if sol_fama.problem
       error('Could not solve the problem with FAMA - infeasible?')
    end 
end

% solve chambolle-pock I
if strcmp(settings.method, 'CPI_G')
    
    % solve the problem in matlab and prepare data for export
    [sol_CPI, stats_CPI, stats_CPI_plot, data_to_export_CPI] = CPI_G(prob, settings);
    
    % export data to a header file
    export_matrices_CPI(data_to_export_CPI);
        
    if sol_CPI.problem
       error('Could not solve the problem with CPI - infeasible?')
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

% CHAMBOLLE-POCKI
if strcmp(settings.method, 'CPI_G')
    semilogy([stats_CPI_plot.rDual], 'b');
    hold on; grid on
    semilogy([stats_CPI_plot.rPrimal], 'b-.');
    legend('Dual residual Chambolle-Pock I', 'Primal residual Chambolle-Pock I');
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

%% Compare to Yalmip

yal_x = sdpvar(n, N);
yal_u = sdpvar(m, N-1);
yal_P = sdpvar(n, n, 'symmetric');

con = [];
for i = 1:N-1
  con = con + set(yal_x(:,i+1) == A*yal_x(:,i) + B*yal_u(:,i));
end

con = con + set( -5 <= yal_x <= 5) + set(-10 <= yal_u <= 10);
con = con + set( norm(0.5*yal_x(:,2) + 2, 2) + yal_x(1,2)/4 <= 100*yal_x(2,2) );
% con = con + set( norm(yal_u(:,3),Inf) <= 2);
con = con + set(yal_x(:,end)'*yal_x(:,end) <= 0.2);
con = con + set(yal_x(:,1) == x0);

yal_obj = 0;
for i = 1:N-1
 yal_obj = yal_obj + yal_x(:,i)'*Q*yal_x(:,i) + yal_u(:,i)'*R*yal_u(:,i);
end
yal_obj = yal_obj + yal_x(:,end)'*yal_x(:,end);
yal_obj = yal_obj + 10.25*norm(yal_x(:,3), 2);

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
cprintf('_blue', '                OPTIMAL SOLUTION from FAMA               \n')
fprintf('x = \n'); disp(full(x.val))
fprintf('u = \n'); disp(full(u.val))
 
cprintf('_blue', '                ERROR TO YALMIP for FAMA              \n')
fprintf('%15s = %.2e\n', '||x - yal_x||', norm(x.val-double(yal_x)));
fprintf('%15s = %.2e\n', '||u - yal_u||', norm(u.val-double(yal_u)));

end

% chambolle-pockI
if strcmp(settings.method, 'CPI_G') 

splitProb.setSolution(prob, sol_CPI);
% Solution is now stored in the variables
cprintf('_blue', '                OPTIMAL SOLUTION from CHAMBOLLE-POCKI      \n')
fprintf('x = \n'); disp(full(x.val))
fprintf('u = \n'); disp(full(u.val))
 
cprintf('_blue', '                ERROR TO YALMIP for CHAMBOLLE-POCKI        \n')
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