%% Infinite horizon MPC using ADMM - main
clear

clear all;
close all;
% randn('state',0);
% rand('state',0);

splitProb.clearProblem;

MAXITER = 1e4;

%% Defining the problem parameters
n=5;

A = 10*randn(n,n); b = 5*randn(n,1);   c = 2.5*randn(n,1); d = 2;
H = eye(n);

%% solve by yalmip

yal_z = sdpvar(n, 1);
con = [];
con = set(cone(A*yal_z + b, c'*yal_z + d));
% con = set(cone(0.5*yal_z + b, 0.2*yal_z + d));
% con = set(norm(10*yal_z + b, 2) <= 10*yal_z + d);
yal_obj = .5*yal_z'*H*yal_z;

result = solvesdp(con, yal_obj, sdpsettings('verbose', 0, 'solver', 'cplex'));

if result.problem, error('YALMIP could not solve problem'); end

z = splitvar(n,1);

obj = .5*z'*H*z;
norm(A*z + b, 2) <= c'*z + d;
minimize(obj);

prob = splitProb.genProblem;

settings.dualTol   = 1e-4; 
settings.primalTol = 1e-4;
settings.maxItr    = 1e5;
settings.method = 'admm';

[sol_admm, stats_admm, data_to_export_admm] = admm(prob, settings);