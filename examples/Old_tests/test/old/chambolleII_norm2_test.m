% cP illustration
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
n = 6;

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
obj = .5*x'*Q*x + l*norm(T*x+t, 2);
A*x == b;
minimize(obj);

% Generate the problem data
prob = splitProb.genProblem;

K1 = [full(prob.prox.L); prob.dat.A];
k1 = [full(prob.prox.l); prob.dat.b];
sigmaG = min(eig(Q));
Mu = max(svds(K1)); % Bound on K operator
theta = 1; sigma = 1*.9/Mu; tau = (1/1)*sigma;
p1 = zeros(size(prob.prox.L, 1), 1); old.p1 = p1;   tilde.p1 = p1;
nu1 = zeros(size(prob.dat.A, 1), 1); old.nu1 = nu1; tilde.nu1 = nu1;
z1 = zeros(size(prob.dat.Q, 1),  1); old.z1 = z1;   bar.z1 = z1;

% precompute some stuff man
matrix_inverse = inv(Q + eye(size(Q,1))*(1/tau));

% Main
for jjj = 1:maxiter 
    
    % step 0: initialization
    old.p1  = p1;
    old.nu1 = nu1;
    old.z1 = z1;
    old.barz1 = bar.z1;
    sigma_prev = sigma;
    tau_prev = tau;
    
    % step 1: compute vector we are going to project
    tilde.p1  = p1 + sigma_prev * (prob.prox.L * old.barz1 + prob.prox.l);
    nu1 = nu1 + sigma_prev * (prob.dat.A * old.barz1 - prob.dat.b);

    p1 = proj_normBall(tilde.p1, prob.prox.dat.p, prob.prox.weight);
    
    % step 2
    matrix_inverse = inv(Q + eye(size(Q,1))*(1/tau_prev));
    z1 = matrix_inverse * ((1/tau)*z1-(K1'*[p1;nu1]));

    % theta, tau, sigma update
    theta = 1/sqrt(1+2 * sigmaG * tau_prev);
    tau = theta * tau_prev;
    sigma = sigma_prev/theta;
    
    % step 3
    bar.z1 = z1 + theta*(z1-old.z1);

    res.s1  = 1/tau*theta * ( z1-old.z1);
    res.r1  = theta*K1*(bar.z1-old.z1) + 1/sigma * ([p1;nu1]-[old.p1; old.nu1]);
    rDual(jjj) = norm(res.s1);    rPrimal(jjj) = norm(res.r1);
    
    if (rDual(jjj) <= 1e-4 && rPrimal(jjj) <= 1e-4)
        break;
    end
    
end

%--------CVX---------%%
cvx_begin
cvx_solver 'sedumi'
    variables Z(2*n,1);
    obj = .5*Z'*Q*Z + l*norm(T*Z+t, 2);
    minimize ( obj )
    subject to 
      A*Z == b;
cvx_end       