clear all;  close all;

N_carts = 2;
%continious ss system
sys = model_generator(N_carts);


%%
Ts = 0.05;

sys_dis = c2d(sys, Ts);
 
A = sys_dis.A; 
B = sys_dis.B;
 C = sys_dis.C;
 D = sys_dis.D;

nx = size(A,1);
nu = size(B,2);

umin = -0.7;
umax = 0.7;


%% SPLIT - parameter settings
kkk = 1; %:length(fcnList_MPC) { @PrFAMA, @AdPrADMM, @AdPrFADMM, @FPDA};
ttt =  7; %horizon N 30 

Q = eye(nx);
R = 0.1*eye(nu);

z0 = 0.5*ones(nx,1);
   
umin = -0.7;
umax = 0.7;



%% algorithm
    
    fcnList_MPC = { @PrFAMA, @AdPrADMM, @AdPrFADMM, @FPDA};
    
 
%% Solve with FoMs
% Set the problem in split format

% Feasibility check
%             if ~(strcmp(cvx_status,'Solved') || strcmp(cvx_status,'Inaccurate/Solved'))
%                 ITER_COUNT(:,countVars,countHorzs) = 0;
%                 break;
%             end

splitProb.clearProblem;
x = splitvar(nx,ttt+1); % declare states as optimization variables
u = splitvar(nu,ttt);   % declare input as optimization variables
x_par = parameter(nx,1); % declare initial state as a paramter

obj = 0; % set objective to zero

% quadratic cost up to horizon N
for j = 1:ttt
    obj = obj + 0.5*x(:,j)'*Q*x(:,j)+0.5*u(:,j)'*R*u(:,j); %   + alpha{iii,ttt}(1)*norm(u(:,j),1)+alpha{iii,ttt}(2)*norm(x(:,j),2);
end

% last cost 
obj = obj + 0.5*x(:,ttt+1)'*Q*x(:,ttt+1);  %+alpha{iii,ttt}(2)*norm(x(:,ttt+1),2);

% assign initial parameter to first time step
x(:,1) == x_par;

% Define constraints
for j = 1:ttt
    x(:,j+1) == A*x(:,j) + B*u(:,j);
    u(:,j) <= umax;
    u(:,j) >= umin;
    %                 norm(FF*x(:,j+1),2) <= FF_x_hs;
    %norm(FF_u*u(:,j),2) <= FF_u_hs;
    %                norm(x(:,j+1),2) <= 2;
    %               norm(u(:,j),2) <= 2;
end

% minimize objective
minimize(obj);

prob = splitProb.genProblem;
settings.maxItr  = 50;
settings.dualTol = 1e-2;
settings.primalTol = 1e-2;


x_par.set(z0); % initial state some random value.

%%% settings
if (kkk==1) % FAMA
    %settings.restart = 'yes';
    %     settings = rmfield(settings, 'restart');
    settings.precond = 'no';
    settings.Mat_Vec = 'FPGA_matvec'; % FPGA_matvec ss
    settings.Lin_Solve = 'auto_FPGA'; % ldl_ss, invert, ldl, ldl_lp, invert_FPGA
    settings.FPGA_PL = 1;
    sd = coder_ama(prob,settings);
    
elseif (kkk==2) % ADMM
    settings.relaxation  = 1.8;
    settings.rho = 1;
    %         settings = rmfield(settings, 'adapt');
    settings.Mat_Vec = 'ss';
    settings.Lin_Solve = 'ldl_lp'; % ldl_ss, invert, ldl, ldl_lp
    sd = coder_admm(prob,settings);
    
elseif (kkk==3) % FADMM
    settings.rho = 1;
    settings.Mat_Vec = 'ss';
    settings.Lin_Solve = 'ldl_lp'; % ldl_ss, invert, ldl, ldl_lp
    sd = coder_admm(prob,settings);
    
elseif (kkk==4) % FPDA
    %settings = rmfield(settings, 'rho');
    settings.Mat_Vec = 'ss';
    settings.Lin_Solve = 'ldl_lp'; % 'auto', 'invert', 'llt', 'ldl_ss', 'ldl_lp', 'llt_ss'
    sd = coder_FPDA(prob,settings);
end


% Solve
[sol, stats, stats_plot] = fcnList_MPC{kkk}(prob, settings);

fom.x = reshape(sol.x(1:(ttt+1)*nx,1),nx,ttt+1);
fom.u = reshape(sol.x(nx*(ttt+1)+1:nx*(ttt+1)+ttt*nu,1),nu,ttt);

sd.add_var('par_ex',z0);
sd.add_var('sol_x',  sol.x);
sd.add_var('sol_lam',  sol.lam);

if (kkk ~=4)
    sd.add_var('sol_y',  sol.y);
end

sd.write_to_file('user_probData');
%             write_FPGA(sd);
% system('make clean; make')
% system('./FAMA')
%system('valgrind --tool=callgrind ./fadmm')


% Stats
% fom.x{kkk,countVars,countHorzs} = reshape(sol.x(1:(ttt+1)*nx,1),nx,ttt+1);
% fom.u{kkk,countVars,countHorzs} = reshape(sol.x(nx*(ttt+1)+1:nx*(ttt+1)+ttt*nu,1),nu,ttt);
% ITER_COUNT(kkk,countVars,countHorzs) = stats.numiter;
% STATS_PLOT{kkk,countVars,countHorzs} = stats_plot;
% ERROR(kkk,countVars,countHorzs) = (norm([fom.u{kkk,countVars,countHorzs}] - U) + norm([fom.x{kkk,countVars,countHorzs}] - X))/(norm(U)+norm(X));
% %             ITERS(kkk,iii,ttt) = mean(ITER_COUNT(kkk,:,:));
% 
% countHorzs = countHorzs+1;
% 
% countVars = countVars+1;


%% FPGA code generation

nPrimal = size(prob.coderData.A,2);
nDual = size(prob.coderData.L,1);
nParam = size(prob.coderData.pL,2);

settings_fpga.algorithm = 'fama';
settings.FPGA_PL = 1;
settings_fpga.x0 = z0;
 settings.resi_tol = [1e-2, 1e-2, 50 , 1];

PL_template_generator(nParam, nPrimal, nDual, settings_fpga);
SPLIT_template

