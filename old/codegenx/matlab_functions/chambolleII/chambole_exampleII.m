clear

clear all;
close all;

% Must always call clearProblem before starting! 
% (or use *clear all*, clear is not enough)
splitProb.clearProblem;

N = 5;
n = 4;
% p = 2;
m = 2;
% [A,B,C,D] = ssdata(drss(n,p,m));

A  = [0.7 -0.1 0.0 0.0; 0.2 -0.5 0.1 0.0; 0.0 0.1  0.1 0.0; 0.5 0.0  0.5 0.5];
B  = [0.0 0.1; 0.1 1.0; 0.1 0.0; 0.0 0.0];
x0 = [ 0.0230; 2.5925; 3.2206; 1.9639 ];

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

-5 <= x <= 5;
-1 <= u <= 1;

% SOC constraint
% norm(0.5*x(:,2) + 2, 2) + x(1,2)/4 <= 10*x(2,2);

% Ellipsoidal constraint
% x(:,end)'*x(:,end) <= 5;
% 0.5*x(:,2) + 2
% 1-norm constraint
% norm(u(:,3),1)*5 + x(:,3)'*4*x(:,3) <= 10;

% x0 = 1.5*ones(n,1);
x(:,1) == x0;

minimize(obj);
% minimize(obj + 3*norm(x(:,3), 1));

%% Generate the problem data
prob = splitProb.genProblem;

%% Solve the problem

settings.dualTol   = 1e-5; 
settings.primalTol = 1e-5;
settings.maxItr    = 1e4;

% NOW I WORK on chambolleI - I think it is ready
% possible values: admm, ama, fadmm, fama 
% settings.method = 'fama';

dat   = prob.dat;
prox  = prob.prox;
nProx = length(prox);

% Validate inputs
% validate_inputs(prob, settings);

% -----------------------
% Pre-processing
% -----------------------

% n = size(dat.Q,1);
% m = size(dat.A,1);

% Setup the KKT matrix and pre-solve
L = []; l = []; 
nProxVars = zeros(nProx,1);
kk = 1;
for i = 1:nProx
  l            = [l;prox(i).l];
  L            = [L;prox(i).L];
  proxInd{i}   = kk + [0:size(prox(i).L,1)-1];
  kk           = kk + size(prox(i).L,1);
  nProxVars(i) = size(prox(i).L,1);
end

% memory allocation for the newly computed variables
p    = zeros(size(prox.L, 1) ,1);
v    = zeros(size(dat.A, 1)  ,1);
z    = zeros(size(prox.L, 2) ,1);
zbar = zeros(size(prox.L, 2) ,1);

% compute the norm of matrix A
normA = max(svds(dat.A));

% this parameter is required to update theta
sigmaG = min(eig(dat.Q));

% two parameters required for convergence
tau = 0.5;
sig = 0.5;

convergence_param = tau * sig * normA;
if(convergence_param < 1)
    fprintf('Convergence can be guaranteed')
end

% chambolles method
for k = 1:settings.maxItr 
        
    % store previous values
    p_prev = p;
    v_prev = v;
    z_prev = z;
    zbar_prev = zbar;
    sig_prev = sig;
    tau_prev = tau;
    
    % compute p according to the formula and compute the projection
    p = p + sig_prev * (prox.L * zbar_prev + prox.l);
    p = min(0, p);
       
    % update the dual variable associated with the dynamic of the system
    v = v + sig_prev * (dat.A * zbar_prev - dat.b);

    % z - update
    temp_matrix = inv(dat.Q + 1/tau_prev * eye(size(dat.Q)));
    z = temp_matrix * (1/tau_prev * z_prev + dat.c' - prox.L' * p - dat.A' * v);
    
    % theta, tau, sigma update
    theta = 1/sqrt(1+2 * sigmaG * tau_prev);
    tau = theta * tau_prev;
    sig = sig_prev/theta;
    
    % zbar update
    zbar = z + theta*(z - z_prev);
    
end

yal_x = sdpvar(n, N);
yal_u = sdpvar(m, N-1);
% yal_P = sdpvar(n, n, 'symmetric');

con = [];
for i = 1:N-1
  con = con + set(yal_x(:,i+1) == A*yal_x(:,i) + B*yal_u(:,i));
end

% con = con + set( -5 <= yal_x <= 5) + set(-1 <= yal_u <= 1);
% con = con + set(yal_x(:,end)'*yal_x(:,end) <= 5);
con = con + set(norm(0.5*x(:,2) + 2, 2) + x(1,2)/4 <= 10*x(2,2));
con = con + set(yal_x(:,1) == x0);

yal_obj = 0;
for i = 1:N-1
 yal_obj = yal_obj + yal_x(:,i)'*Q*yal_x(:,i) + yal_u(:,i)'*R*yal_u(:,i);
end
yal_obj = yal_obj + yal_x(:,end)'*yal_x(:,end);
% yal_obj = yal_obj + 3*norm(yal_x(:,3), 1);

result = solvesdp(con, yal_obj, sdpsettings('verbose', 0, 'solver', 'cplex'));

if result.problem, error('YALMIP could not solve problem'); end

% reshape yalmip solution for comparison and create one array
yal_x_reshaped = reshape(double(yal_x), n*N, 1);
yal_u_reshaped = reshape(double(yal_u), m*(N-1), 1);
yal_solution   = [ yal_x_reshaped; yal_u_reshaped ];

% error yalmip vs matlab
err = norm(z - yal_solution);

fprintf('\nError yalmip vs matlab: %.2g\n', err);