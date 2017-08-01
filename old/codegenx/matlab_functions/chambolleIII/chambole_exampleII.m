problemData.n = 4;
problemData.m = 2;
problemData.N = 4;
% problemData.c = ones(problemData.n + problemData.m, 1);
problemData.Q = eye(problemData.n);
problemData.R = eye(problemData.m);

% system dynamics
dynamics.A = [0.7 -0.1 0.0 0.0; 0.2 -0.5 0.1 0.0; 0.0 0.1  0.1 0.0; 0.5 0.0  0.5 0.5];
dynamics.B = [0.0 0.1; 0.1 1.0; 0.1 0.0; 0.0 0.0];
dynamics.x0 = rand(problemData.n, 1)*5.25;

% dynamics.A  = rand(problemData.n, problemData.n);
% dynamics.B  = rand(problemData.n, problemData.m);
% dynamics.x0 = rand(problemData.n, 1)*2.25;

% polytopic constraints
constraint.C = [eye(problemData.n); -eye(problemData.n)];
constraint.D = [eye(problemData.m); -eye(problemData.m)];
constraint.c = ones(2*problemData.n, 1)*10;
constraint.d = ones(2*problemData.m, 1);
constraint.type = 'box';

% obtain solution 
[solver_solution, solution] = obtain_solution(dynamics, constraint, problemData);

% create example for chambolle - modify the file later
[probSetup, proxSetup] = create_examples_for_chambolle(dynamics, constraint, problemData);

% small modification for sake of simlicity
proxSetup.EN = proxSetup.EN';
probSetup.c  = probSetup.c';
probSetup.c = reshape(probSetup.c, 1, 4*6);

% declare parameters specific for Chambolle's method
% memory allocation and initilization
p_prev = zeros(size(proxSetup.en,1), 1);
v_prev = zeros(size(probSetup.b,1), 1);
zbar_prev = zeros(size(proxSetup.EN,2), 1);
z_prev = zbar_prev;


% memory allocation for the newly computed variables
p = zeros(size(p_prev));
v = zeros(size(v_prev));
z = zeros(size(z_prev));
zbar = zeros(size(zbar_prev));

% compute the norm of matrix A
normA = max(svds(probSetup.A));

% this parameter is required to update theta
sigmaG = min(eig(probSetup.Q));

% two parameters required for convergence
tau_prev = 0.5;
sig_prev = 0.5;

convergence_param = tau_prev * sig_prev * normA;
if(convergence_param < 1)
    fprintf('Convergence can be guaranteed')
end

% variable for the maximal number of iteration
maxiter = 1e4;

% chambolles method
for k = 1:maxiter
    
    % compute p according to the formula and compute the projection
    p = p + sig_prev * (proxSetup.EN * zbar_prev - proxSetup.en);
    p = max(0, p);
       
    % update the dual variable associated with the dynamic of the system
    v = v + sig_prev * (probSetup.A * zbar_prev - probSetup.b);

    % z - update
    temp_matrix = inv(probSetup.Q + 1/tau_prev * eye(size(probSetup.Q)));
    z = temp_matrix * (1/tau_prev * z_prev + probSetup.c' - proxSetup.EN' * p - probSetup.A' * v);
    
    % theta, tau, sigma update
    theta = 1/sqrt(1+2 * sigmaG * tau_prev);
    tau = theta * tau_prev;
    sig = sig_prev/theta;
    
    % zbar update
    zbar = z + theta*(z - z_prev);
    
    % store previous values
    p_prev = p;
    v_prev = v;
    z_prev = z;
    zbar_prev = zbar;
    sig_prev = sig;
    tau_prev = tau;
    
end












