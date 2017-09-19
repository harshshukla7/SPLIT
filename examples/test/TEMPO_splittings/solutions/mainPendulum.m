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
dat.L = [zeros(2*N*nu,(N+1)*nx) kron(eye(N), [eye(nu); -eye(nu)])];
dat.l = kron(ones(N,1),[umax;-umin]);
L = dat.L;
l = dat.l;


SimTime = Tf/Ts;

%% ----------------------------------------------
%% Solve w/ splitting methods - MPC for regulation 
%% ----------------------------------------------
fcnList_MPC = { @AMA, @FAMA, @PrFAMA };

clear settings;
clear ITER_COUNT;

%% Iterate over the algorithms
for kkk = 1:length(fcnList_MPC)
    % Global settings for solver
    settings.maxItr = 150;
    settings.dualTol = 1e-3;
    settings.primalTol = 1e-3;
    if (kkk==1) % AMA w/o preconditioning
    elseif (kkk==2) % FAMA w/o preconditioning
        settings.restart = 'yes';
    elseif (kkk==3) % FAMA w/ preconditioning
        [settings.dualPrec] = precAMA(dat); 
        settings.rho = 1;
        settings.restart = 'yes';
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
    %% Receding horizon control
    for ttt = 1:SimTime-N-1   
        % Solve
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % FILL IN THE FUNCTIONS
        [sol, stats, stats_plot] =...
            fcnList_MPC{kkk}(dat,dxs,dus,state_sim(:,end),settings);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        x = splitvar(nx,N+1);
        x0 = parameter(nx,1);
        u = splitvar(nu,N);
        
        
        
        % Recover optimal closed loop trajectories
        fom.x{kkk} = reshape(sol.x(1:(N+1)*nx),nx,N+1);
        fama.x{kkk} = [fama.x{kkk} fom.x{kkk}(:,2)];
        fom.u{kkk} = reshape(sol.x((N+1)*nx+1:end),nu,N);
        fama.u{kkk} = [fama.u{kkk} fom.u{kkk}(:,1)];
        
        if (ttt==SimTime-N-1)
            fom.x{kkk} = reshape(sol.x(1:(N+1)*nx),nx,N+1);
            fama.x{kkk} = [fama.x{kkk} fom.x{kkk}(:,2:end)];
            fom.u{kkk} = reshape(sol.x((N+1)*nx+1:end),nu,N);
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
        ITER_COUNT(kkk,ttt) = stats.numiter;
    end
    %% Keep statistics for each algorithm
    clf;
    pause;
    ITERS(kkk) = mean(ITER_COUNT(kkk,:));
end