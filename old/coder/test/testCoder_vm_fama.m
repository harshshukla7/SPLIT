clear

splitProb.clearProblem;

x = splitvar(3,1);
z = parameter(2,1);

Q = eye(3);
F = randn(3,2);
K = randn(20,3);
k = ones(20,1);

K*x <= k;
minimize(z'*F'*x + x'*Q*x);

prob = splitProb.genProblem;

% Value of the parameter for testing
par = [1;1];

%% Solve using CPLEX

y = sdpvar(3,1);

ddd = solvesdp(set(K*y <= k), par'*F'*y + y'*Q*y)
y = double(y);

%% Solve using Giorgos' solver




%% Compile a VM FAMA solver
cdr = splitCoder(prob, 'vm_fama','funcVMFAMA');

% Diagonal preconditioning for acquiring VM
cvx_begin sdp
  cvx_solver sedumi
  variable invVM(size(K,1), size(K,1)) diagonal
  minimize ( trace(invVM) )
  subject to
  invVM >= K*(2*eye(3)\K');
cvx_end
VM = inv(full(invVM));


cdr.VM = VM;
cdr.SP   = eye(size(cdr.L',1));
cdr.dualTol = 1e-6;
cdr.primalTol = 1e-6;

cdr.rho = 1;
cdr.itr_conv = 10;
cdr.KKTSolveMode = 'Invert KKT';
cdr.MatVec_SparseLimit = 0.1;

% cdr.genMex
solverVMFAMA = cdr.genMatlab;

%% Compile an ADMM solver

cdr = splitCoder(prob, 'admm','funcADMM');

cdr.rho = 1;
cdr.itr_conv = 10;
cdr.KKTSolveMode = 'Invert KKT';
cdr.MatVec_SparseLimit = 0.1;
cdr.dualTol = 1e-6;
cdr.primalTol = 1e-6;

% cdr.genMex
solverADMM = cdr.genMatlab;


%% Test if it works

yx = sdpvar(3,1);
yy = sdpvar(2,1);

par = 1*randn(2,1);

ddd = solvesdp(set(K*yx <= k), par'*F'*yx + yx'*Q*yx)
yxVal = double(yx);

% funcs = { @funcVMFAMA_mex, solverVMFAMA, @funcADMM_mex, solverADMM };
funcs = { solverVMFAMA, solverADMM };

for i = 1:length(funcs)
  [x_(:,i), y_(:,i), lambda_(:,i), itr(i), tm(i)] = funcs{i}(par);
end

%%
splitProb.clearProblem;
x = splitvar(3,1);
K*x <= k;
minimize(par'*F'*x + x'*Q*x);

prob = splitProb.genProblem;
settings.variable_metric = VM;
[sol, stats, stats_plot, ~] = vm_fama_MPC(prob, settings);
xx = sol.x;


K*yxVal - k

[x_ yxVal]