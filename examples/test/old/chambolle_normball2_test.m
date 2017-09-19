clear
clear all
close all

% fix the seed
randn('state',0);
rand('state',0);

% script to test chambolle

% clear everything
splitProb.clearProblem;

% maximum number of iterations
maxiter = 1e5;

% number of states
n = 5;

% l2 - ball
T = randn(2*n); t = randn(2*n, 1); l = 5.5;

% equality constraint
A = randn(n, 2*n); b = randn(n, 1);

% stuff for the objective
Q = randn(2*n, 2*n);
Q = Q'*Q;

% one splitvar for the problem
x = splitvar(2*n, 1);

% formulate the problem
obj = .5*x'*Q*x;
norm(T*x+t, Inf) <= l;
A*x == b;
minimize(obj);

% Generate the problem data
prob = splitProb.genProblem;

% second formulation - to support the output from the interface
% prepare matrices for the computation
K = [full(prob.prox.L); prob.dat.A];
k = [full(prob.prox.l); prob.dat.b];
Mu = max(svds(K)); % Bound on K operator
theta = 1; sigma = 1*.9/Mu; tau = (1/1)*sigma;
p1 = zeros(size(prob.prox.L, 1), 1); old.p1 = p1;   tilde.p1 = p1;
nu1 = zeros(size(prob.dat.A, 1), 1); old.nu1 = nu1; tilde.nu1 = nu1;
z1 = zeros(size(prob.dat.Q, 1),  1); old.z1 = z1;   bar.z1 = z1;

% construct everything manually
% K = [T; A];
% k = [t; b];
% Mu = max(svds(K)); % Bound on K operator
% theta = 1; sigma = 1*.9/Mu; tau = (1/1)*sigma;
% p1 = zeros(2*n,1); old.p1 = p1; tilde.p1 = p1; % has the size of T
% nu1 = zeros(n,1); old.nu1 = nu1; tilde.nu1 = nu1;
% z1 = zeros(2*n,1);  old.z1 = z1;  bar.z1 = z1;

% precompute some stuff man
matrix_inverse = inv(Q + eye(size(Q,1))*(1/tau));

for jjj = 1:maxiter
    
    % step 0: initialize the stuff
    old.p1 = p1;
    old.nu1 = nu1;
    
    % step 1: compute vector we are going to project
    tilde.p1  = p1 + sigma*(prob.prox.L * bar.z1 + prob.prox.l);
    nu1 = nu1 + sigma*(prob.dat.A * bar.z1 - prob.dat.b);
    
    % compute projection
    p1 = prox_norm(tilde.p1, prob.prox.dat.c*sigma, prob.prox.dat.pDual);
    
    % step 2
    old.z1 = z1;
    z1 = matrix_inverse * ((1/tau)*z1-(K'*[p1;nu1]));

    % step 3
    bar.z1 = z1 + theta*(z1-old.z1);
    
    res.s1  = 1/tau*theta*(z1-old.z1);
    res.r1  = theta*K*(bar.z1-old.z1) + 1/sigma*([p1;nu1]-[old.p1;old.nu1]);
            
    rDual(jjj) = norm(res.s1);    rPrimal(jjj) = norm(res.r1);
    if (rDual(jjj) <= 1e-4 && rPrimal(jjj) <= 1e-4)
        break;
    end
        
end

% solve the problem by cvx
cvx_begin
cvx_solver 'sedumi'
    variables Z(2*n,1);
    obj = .5*Z'*Q*Z;
    minimize ( obj )
    subject to 
      A*Z == b;
      norm(T*Z+t, Inf) <= l;
cvx_end