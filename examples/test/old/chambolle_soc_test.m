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

% soc constraint
T = randn(3*n,2*n);   t = .18 - randn(3*n,1);
c = randn(2*n,1); d = 3.9;

% equality constraint
A = randn(n,2*n); b = randn(n,1);

% stuff for the objective
Q = randn(2*n, 2*n);
Q = Q'*Q;

% one splitvar for the problem
x = splitvar(2*n, 1);

% formulate the problem
obj = .5*x'*Q*x;
norm(T*x + t, 2) <= c'*x + d;
A*x == b;
minimize(obj);

% Generate the problem data
prob = splitProb.genProblem;

%% first formulation

% % construct everything manually
% K = [[T;c']; A];
% k = [[t;d]; b];
% no.con.all = size(K,1);
% Mu = max(svds(K)); % Bound on K operator
% theta = 1; sigma = 1*.9/Mu; tau = (1/1)*sigma;
% p = zeros(3*n+1,1); old.p = p; tilde.p = p; % has the size of T
% nu = zeros(n,1); old.nu = nu; tilde.nu = nu;
% z = zeros(2*n,1);   old.z = z;  bar.z = z;
% 
% % precompute some stuff man
% matrix_inverse = inv(Q + eye(size(Q,1))*(1/tau));
% 
% for jjj = 1:maxiter
%     
%     % step 0: initialize the stuff
%     old.p = p;
%     old.nu = nu;
%     
%     % step 1: compute vector we are going to project
%     tilde.p = p + [sigma*(T*bar.z+t); sigma*(c'*z+d)];
%     nu = nu + sigma*(A*bar.z-b);
%     
%     % compute projection
%     p = proj_secondOrderConeConj(tilde.p);
%     
%     % step 2
%     old.z = z;
%     z = matrix_inverse * ((1/tau)*z-(K'*[p;nu]));
% 
%     % step 3
%     bar.z = z + theta*(z-old.z);
% 
%     res.s  = 1/tau*theta*(z-old.z);
%     res.r  = theta*K*(bar.z-old.z) + 1/sigma*([p;nu]-[old.p;old.nu]);
%             
%     rDual(jjj) = norm(res.s);    rPrimal(jjj) = norm(res.r);
%     if (rDual(jjj) <= 1e-4 && rPrimal(jjj) <= 1e-4)
%         break;
%     end
% end

%% second formulation - to support the output from the interface

% construct everything manually
K = [[c';T]; A];
k = [[d;t]; b];
Mu = max(svds(K)); % Bound on K operator
theta = 1; sigma = 1*.9/Mu; tau = (1/1)*sigma;
p1 = zeros(3*n+1,1); old.p1 = p1; tilde.p1 = p1; % has the size of T
nu1 = zeros(n,1); old.nu1 = nu1; tilde.nu1 = nu1;
z1 = zeros(2*n,1);  old.z1 = z1;  bar.z1 = z1;

% precompute some stuff man
matrix_inverse = inv(Q + eye(size(Q,1))*(1/tau));

for jjj = 1:maxiter
    
    % step 0: initialize the stuff
    old.p1 = p1;
    old.nu1 = nu1;
    
    % step 1: compute vector we are going to project
    tilde.p1  = p1 + [sigma*(c'*z1+d); sigma*(T*bar.z1+t)];
    nu1 = nu1 + sigma*(A*bar.z1-b);
    
    % compute projection
    p1 = proj_secondOrderConeConj(tilde.p1);
    
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
      { T*Z + t,c'*Z + d } == lorentz(3*n);
cvx_end