
clear all;
close all;
clc;

XINIT = [0 0.12 0 0].';
Ts = 0.05;
Tf = 7;
input.Ts = Ts;
input.nSteps = 3;

% load data
load datPendulum
nx = dat.dim.nx;
nu = dat.dim.nu;
N = dat.time.N;
dxs = zeros((N+1)*nx,1);
dus = zeros(N*nu,1);
Q = dat.cost.Q;
R = dat.cost.R;
P = dat.cost.P;
A = dat.model.A;
B = dat.model.B;

% constraints
umax = 0.75;
umin = -3.25;

splitProb.clearProblem;
% Write tracking optimization problem in Split format
x = splitvar(nx,N+1);
u = splitvar(nu,N);
% z=(x,u)
zs = parameter(nx*(N+1)+nu*N+nx,1);
obj = 0;
for jjj = 1:N
    obj = obj + (x(:,jjj)-zs((jjj-1)*nx+1:jjj*nx))'*Q*(x(:,jjj)-zs((jjj-1)*nx+1:jjj*nx)) +...
        (u(:,jjj)-zs((N+1)*nx+(jjj-1)*nu+1:(N+1)*nx+jjj*nu))'*R*(u(:,jjj)-zs((N+1)*nx+(jjj-1)*nu+1:(N+1)*nx+jjj*nu));
    u(:,jjj) <= umax;
    u(:,jjj) >= umin;
end
obj = obj + (x(:,N+1)-zs(N*nx+1:(N+1)*nx))'*P*(x(:,N+1)-zs(N*nx+1:(N+1)*nx));
x(:,1) == zs(end-nx+1:end,1);
x(:,2:end) == A*x(:,1:end-1) + B*u;
minimize(obj);
% Generate the problem data
prob = splitProb.genProblem;

SimTime = Tf/Ts;

%% mex setup :
opt.primalTol = 1e-3;
opt.dualTol   = 1e-3;
opt.MAXITR    = 150;
opt.ITR_PER_CONV_TEST = 5;

%If you want to time just the solve-time, then set this larger than one,
%and the time returned will be the average per solve
opt.number_of_solves = 1;

%%
fcnList_MPC = { @PrFAMA, @PrFAMA, @PrFAMA};
clear settings;
clear ITER_COUNT;

%% Iterate over the algorithms

for kkk = 3 % :length(fcnList_MPC)
    % Global settings for solver
    settings.maxItr = 150;
    settings.dualTol = 1e-3;
    settings.primalTol = 1e-3;
    if (kkk==1) % AMA w/o preconditioning
        %% Generate for data for c
        clear tp_settings
        tp_settings.precond = 'no';
        tp_settings.Mat_Vec = 'ss';
        tp_settings.Lin_Solve = 'ldl_ss';
        clear sd_ama;
        sd_ama = coder_ama(prob,tp_settings);
    elseif (kkk==2) % FAMA w/o preconditioning
        
        settings.restart = 'yes';
        
        %% Generate for data for c
        clear tp_setting
        tp_settings.precond = 'no';
        tp_settings.Mat_Vec = 'ss';
        tp_settings.Lin_Solve = 'ldl_ss';
        tp_settings.adaptive_restart = 'yes';
        clear sd_ama;
        sd_ama = coder_ama(prob,tp_settings);
        
    elseif (kkk==3) % FAMA w/ preconditioning
        [settings.scale_rDual, ~] = equilPrecDual_ADMM_AMA(prob, 'fama', 'equilibration');
        settings.rho = 1;
        settings.restart = 'yes';
        
        %% Generate for data for c
        clear tp_setting
        tp_settings.precond = 'yes';
        tp_settings.Mat_Vec = 'ss';
        tp_settings.Lin_Solve = 'ldl_ss';
        tp_settings.E = (settings.scale_rDual);
        tp_settings.D_scale = eye(size(prob.dat.Q,1));
        tp_settings.rho = 1;
        tp_settings.adaptive_restart = 'yes';
        
        clear sd_ama;
        sd_ama = coder_ama(prob,tp_settings);
    end
    
    % save optimal x and u for each method
    fom.x{kkk} = [];
    fama.x{kkk,1} = XINIT;
    fom.u{kkk} = [];
    fama.u{kkk} = [];
    
    state_sim = XINIT;
    state_lin = XINIT;
    iter = 0;
    time = 0;
    
    
    
    %% generate split data
    
%     clear sd_ama;
%     sd_ama = coder_ama(prob,tp_settings);
    
    sd_ama.write_to_file('probData');
    
    
    %% use mex file to compile
    
    mex -output ama_mex  -v  -largeArrayDims  fama.c matrix_ops.c probData.c  splitLoad.c splitTimer.c ama_mex.c ldl.c


    %% Receding horizon control
    
    
    for ttt = 1 :SimTime-N-1
        
        % set the parameter
        %     clear zs.set
        
        par_ex = [dxs;dus;state_sim(:,end)];
        zs.set(par_ex);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % FILL IN THE FUNCTIONS
        %[sol, stats, stats_plot] = fcnList_MPC{kkk}(prob, settings);
        
        %% write and similate in c
        
        sol = ama_mex(par_ex, opt);
        
        % Recover optimal closed loop trajectories
        fom.x{kkk} = reshape(sol.primal(1:(N+1)*nx),nx,N+1);
        fama.x{kkk} = [fama.x{kkk} fom.x{kkk}(:,2)];
        fom.u{kkk} = reshape(sol.primal((N+1)*nx+1:end),nu,N);
        fama.u{kkk} = [fama.u{kkk} fom.u{kkk}(:,1)];
        
        if (ttt==SimTime-N-1)
            fom.x{kkk} = reshape(sol.primal(1:(N+1)*nx),nx,N+1);
            fama.x{kkk} = [fama.x{kkk} fom.x{kkk}(:,2:end)];
            fom.u{kkk} = reshape(sol.primal((N+1)*nx+1:end),nu,N);
            fama.u{kkk} = [fama.u{kkk} fom.u{kkk}];
        end
        
        % apply control to nonlinear system
        input.x = state_sim(:,end);
        input.u = max(min(fama.u{kkk}(ttt),umax),umin);
        output = RK4_integrator( @ode, input );
        state_sim(:,end+1) = output.value;
        
        % linear system
        state_lin(:,end+1) = A*state_lin(:,end) + B*max(min(fama.u{kkk}(ttt),umax),umin);
        
        % next time step and visualize result
        iter = iter+1;
        time(end+1) = iter*Ts;
        visualize(time, state_sim, 0.8);
        
        pause(Ts/2);
%        ITER_COUNT(kkk,ttt) = stats.numiter;
        
        
        %% compare the error
       %% error_x =  norm(sol_mex.primal-sol.x)


        
    end
    
    %% Keep statistics for each algorithm
    clf;
%    pause;
%    ITERS(kkk) = mean(ITER_COUNT(kkk,:));
    
end


%% mex interface

% mex  -v  -largeArrayDims ama_mex.c fama.c matrix_ops.c probData.c splitLoad.c splitTimer.c TEMPO_testCoder.c  -L/usr/local/lib ...
%     LIBS='-lm -lldl' -lmwlapack -lmwblas -I/usr/local/include LDFLAGS='-L/usr/local/lib'

% mex -output ama_mex  -v  -largeArrayDims  fama.c matrix_ops.c probData.c  splitLoad.c splitTimer.c ama_mex.c 
% 
% opt.primalTol = 1e-3;
% opt.dualTol   = 1e-3;
% opt.MAXITR    = 150;
% opt.ITR_PER_CONV_TEST = 5;
% 
% %If you want to time just the solve-time, then set this larger than one,
% %and the time returned will be the average per solve
% opt.number_of_solves = 1;
% 
% 
% sol_mex = ama_mex(par_ex, opt);

%error_x =  norm(sol_mex.x-sol.x);
%error_x =  norm(sol_mex.dual-sol.lam);
