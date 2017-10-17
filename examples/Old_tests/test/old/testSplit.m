clear

clear all;
close all;
randn('state',1);
rand('state',0);

% Must always call clearProblem before starting! 
% (or use *clear all*, clear is not enough)
splitProb.clearProblem;

N = 5;
n = 3;
p = 2;
m = 2;
[A,B,C,D] = ssdata(drss(n,p,m));

x = splitvar(n, N);
u = splitvar(m, N-1);
x0 = parameter(n,1);
Q = randn(n); Q = Q*Q';
R = randn(m); R = R*R';

obj = 0;
for i = 1:N-1
  x(:,i+1) == A*x(:,i) + B*u(:,i);
  obj = obj + x(:,i)'*Q*x(:,i) + u(:,i)'*R*u(:,i);
end
obj = obj + x(:,end)'*x(:,end);
obj = obj + 5.5*norm(x(:,2), 2);

-5  <= x <= 5;
-10 <= u <= 10;

% % SOC constraint
norm(0.5*x(:,2) + 2, 2) + x(1,2)/4 <= 10*x(2,2);

% ellipsoidal constraint
x(:,end)'*x(:,end) <= 0.0008;

% 1-normball constraint
norm(u(:,1), 2) <= 1.5;

% x0 = 1.5*ones(n,1);
x(:,1) == x0;

minimize(obj);

%% Generate the problem data
prob = splitProb.genProblem;

%% Solve the problem

settings.dualTol   = 1e-5; 
settings.primalTol = 1e-5;
settings.maxItr    = 1e5;
% settings.rho       = 1;

% possible values: admm, ama, fadmm, fama, CPI, CPII
settings.method = 'CPI';

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

% solve CPI
if strcmp(settings.method, 'CPI')
    
    % solve the problem in matlab and prepare data for export
    [sol_CPI, stats_CPI, stats_CPI_plot, data_to_export_CPI] = CPI(prob, settings);
    
    % export data to a header file
    export_matrices_CHAMBOLLEI(data_to_export_chambolle1);
        
    if sol_CPI.problem
       error('Could not solve the problem with CPI - infeasible?')
    end 
end

% solve CPII
if strcmp(settings.method, 'CPII')
    
    % solve the problem in matlab and prepare data for export
    [sol_CPII, stats_CPII, stats_CPII_plot, data_to_export_CPII] = CPII(prob, settings);
    
    % export data to a header file
    % export_matrices_CHAMBOLLEIa(data_to_export_chambolle1);
        
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
    semilogy([stats_ama_plot.rPrimal], 'b-.');
    hold on; grid on
    legend('Dual residual AMA', 'Primal residual AMA');
end

% FAMA
if strcmp(settings.method, 'fama')
    semilogy([stats_fama_plot.rDual], 'b');
    hold on; grid on
    semilogy([stats_fama_plot.rPrimal], 'b-.');
    legend('Dual residual FAMA', 'Primal residual FAMA');
end

% CHAMBOLLE1
if strcmp(settings.method, 'CPI')
    semilogy([stats_CPI_plot.rDual], 'b');
    hold on; grid on
    semilogy([stats_CPI_plot.rPrimal], 'b-.');
    legend('Dual residual Chambolle1', 'Primal residual Chambolle1');
end

% CHAMBOLLE2
if strcmp(settings.method, 'CPII')
    semilogy([stats_CPII_plot.rDual], 'b');
    hold on; grid on
    semilogy([stats_CPII_plot.rPrimal], 'b-.');
    legend('Dual residual Chambolle2', 'Primal residual Chambolle2');
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
con = con + set(norm(0.5*yal_x(:,2) + 2, 2) + yal_x(1,2)/4 <= 10*yal_x(2,2));
con = con + set(norm(yal_u(:,1), 2) <= 1.5);
con = con + set(yal_x(:,end)'*yal_x(:,end) <= 0.0008);
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

% chambolle1
if strcmp(settings.method, 'CPI') 

splitProb.setSolution(prob, sol_chambolle1);
% Solution is now stored in the variables
cprintf('_blue', '                OPTIMAL SOLUTION from CPI      \n')
fprintf('x = \n'); disp(full(x.val))
fprintf('u = \n'); disp(full(u.val))
 
cprintf('_blue', '                ERROR TO YALMIP for CPI        \n')
fprintf('%15s = %.2e\n', '||x - yal_x||', norm(x.val-double(yal_x)));
fprintf('%15s = %.2e\n', '||u - yal_u||', norm(u.val-double(yal_u)));

end

% chambolle2
if strcmp(settings.method, 'CPII') 

splitProb.setSolution(prob, sol_chambolle2);
% Solution is now stored in the variables
cprintf('_blue', '                OPTIMAL SOLUTION from CPII      \n')
fprintf('x = \n'); disp(full(x.val))
fprintf('u = \n'); disp(full(u.val))
 
cprintf('_blue', '                ERROR TO YALMIP for CPII        \n')
fprintf('%15s = %.2e\n', '||x - yal_x||', norm(x.val-double(yal_x)));
fprintf('%15s = %.2e\n', '||u - yal_u||', norm(u.val-double(yal_u)));

end