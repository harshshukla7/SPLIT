clear all;  close all;

N_carts = 4; % number of carts
Ts = 0.5; % sampling time

%continious ss system
sys = model_generator_MSD(N_carts);


%%

sys_dis = c2d(sys, Ts); %discrete ss system

% extract the matrices
A = sys_dis.A;
B = sys_dis.B;
C = sys_dis.C;
D = sys_dis.D;

nx = size(A,1); % number of states
nu = size(B,2); % number of inputs


%% SPLIT - parameter settings

ttt =  10; %horizon N

Q = 50*eye(nx);
R = eye(nu);

z0 = 0.5*ones(nx,1);

umin = -0.5;
umax = 0.5;


%% Soft constraint part

sigma_1 = 1*ones(N_carts,1);

sigma_2 = 1*eye(N_carts);

ri = 0;
x_ci = 0;

xmin = -0.5*ones(N_carts,1);
xmax = 0.5*ones(N_carts,1);


%% Solve with FoMs
% Set the problem in split format


splitProb.clearProblem;
x = splitvar(nx,ttt+1); % declare states as optimization variables
u = splitvar(nu,ttt);   % declare input as optimization variables
x_par = parameter(nx,1); % declare initial state as a paramter
delta = splitvar(N_carts,ttt+1); % for softconstraints

obj = 0; % set objective to zero

% quadratic cost up to horizon N
for j = 1:ttt
    obj = obj + 0.5*x(:,j)'*Q*x(:,j)+0.5*u(:,j)'*R*u(:,j)  + 0.5*delta(:,j)'*sigma_2*delta(:,j) + sigma_1'*delta(:,j);
end

% cost for N+1
obj = obj + 0.5*x(:,ttt+1)'*Q*x(:,ttt+1) + 0.5*delta(:,ttt+1)'*sigma_2*delta(:,ttt+1) ;

% assign initial parameter to first time step
x(:,1) == x_par;


% Define constraints
for j = 1:ttt
    
    x(:,j+1) == A*x(:,j) + B*u(:,j); % dynamics
    u(:,j) <= umax;
    u(:,j) >= umin;
    x(1:2:end,j) >= xmin - delta(:,j);
    x(1:2:end,j) <= xmin + delta(:,j);
    
end

j = ttt+1;
x(1:2:end,j) >= xmin - delta(:,j);
x(1:2:end,j) <= xmin + delta(:,j);

% minimize objective
minimize(obj);

prob = splitProb.genProblem;
settings.maxItr  = 200;
settings.dualTol = 1e-2;
settings.primalTol = 1e-2;



% FAMA
%settings.restart = 'yes';
%     settings = rmfield(settings, 'restart');
settings.precond = 'no';
settings.Mat_Vec = 'ss'; % user sparse matvec
settings.Lin_Solve = 'ldl_ss'; % use sparse ldl factorization
settings.FPGA_PL = 0;
settings.paral = 30;
sd = coder_ama(prob,settings);



x_par.set(z0); % initial state some random value.
sd.add_var('par_ex',z0); % add it for code generation as well (only for processors)

% Solve in matlab
[sol, stats, stats_plot] = PrFAMA(prob, settings);

% extract states and inputs
fom.x = reshape(sol.x(1:(ttt+1)*nx,1),nx,ttt+1);
fom.u = reshape(sol.x(nx*(ttt+1)+1:nx*(ttt+1)+ttt*nu,1),nu,ttt);


% save solution in c file for code generation - (optional)

sd.add_var('sol_x',  sol.x);
sd.add_var('sol_lam',  sol.lam);
sd.add_var('sol_y',  sol.y);


% now write all the data to the file name user_probdata
sd.write_to_file('user_probData');


%% FPGA code generation

nPrimal = size(prob.coderData.A,2);
nDual = size(prob.coderData.L,1);
nParam = size(prob.coderData.pL,2);

settings_fpga.algorithm = 'fama';
settings_fpga.x0 = z0;
settings_fpga.resi_tol = [1e-2, 1e-2, 200 , 200];

Processor_template_generator(nParam, nPrimal, nDual, settings_fpga)
% SPLIT_template

