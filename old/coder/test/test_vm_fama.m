clear

% Simple problem:
%  min x'*Q*x + h'*x
%  s.t. K*x <= k

Q = eye(3);
K = randn(20,3);
k = ones(20,1);
F = 10*randn(2,3);

val = randn(2,1);

%% Solve using CPLEX

y = sdpvar(3,1);
solvesdp(set(K*y <= k), val'*F*y + y'*Q*y)
y = double(y);

%% Solve using Giorgos' solver

splitProb.clearProblem;
x = splitvar(3,1);
K*x <= k;
minimize(val'*F*x + x'*Q*x);
prob = splitProb.genProblem;

% Diagonal preconditioning for acquiring VM
% cvx_begin sdp
%   cvx_solver sedumi
%   variable invVM(size(K,1), size(K,1)) diagonal
%   minimize ( trace(invVM) )
%   subject to
%   invVM >= K*(2*eye(3)\K');
% cvx_end
invVM = sdpvar(size(K,1),1);
solvesdp(diag(invVM) >= K*(2*eye(3)\K'), sum(invVM))
VM = diag(1./double(invVM));
settings.variable_metric = VM;
settings.scale_rPrimal   = eye(size(prob.prox.L',1));
settings.dualTol   = 1e-6;
settings.primalTol = 1e-6;

[sol, stats, stats_plot, ~] = vm_fama_MPC(prob, settings);
xStath = sol.x;

%% Solve using parametric VM_FAMA_MPC

splitProb.clearProblem
x = splitvar(3,1);
par = parameter(2,1);
K*x <= k;
minimize(par'*F*x + x'*Q*x);
prob = splitProb.genProblem;

% Diagonal preconditioning for acquiring VM
% cvx_begin sdp
%   cvx_solver sedumi
%   variable invVM(size(K,1), size(K,1)) diagonal
%   minimize ( trace(invVM) )
%   subject to
%   invVM >= K*(2*eye(3)\K');
% cvx_end
% VM = inv(full(invVM));
invVM = sdpvar(size(K,1),1);
solvesdp(diag(invVM) >= K*(2*eye(3)\K'), sum(invVM))
VM = 1./double(invVM);

cdr = splitCoder(prob, 'vm_fama','funcVMFAMA');

cdr.VM        = VM;
cdr.SP        = ones(size(prob.prox.L',1),1);
cdr.dualTol   = 1e-6;
cdr.primalTol = 1e-6;

cdr.rho = 1;
cdr.itr_conv = 10;
cdr.KKTSolveMode = 'Invert KKT';
cdr.MatVec_SparseLimit = 1;

cdr.genMex
solverVMFAMA = cdr.genMatlab;
[x_, y_, lambda_, itr, tm, rDual] = solverVMFAMA(val, 5, 1e-6);
[x2_, y2_, lambda2_, itr2, tm2, rDual2] = funcVMFAMA_mex(val, 5, 1e-6);


%% Compare solutions

[y xStath x_ x2_]

fprintf('Number of active constraints: %i \n', sum(abs(K*y-k) < 1e-3));

%% Plot of iterations vs dual residual

DUAL = zeros(1000,1);
for i = 1:1000
  [x2_, y2_, lambda2_, itr2, tm2, rDual2] = funcVMFAMA_mex(val, i, 1e-6);
  DUAL(i) = rDual2;
end
