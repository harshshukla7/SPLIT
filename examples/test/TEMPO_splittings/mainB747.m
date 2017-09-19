clear all;  close all;
randn('state', 0);

%% load data and simplify notation
load datB747
nx = dat.dim.nx;
nu = dat.dim.nu;
N = dat.time.N;
SimTime = dat.time.sim;
dxs = dat.ss.x;
dus = dat.ss.u;
umax = dat.constraints.ulims_max;
umin = dat.constraints.ulims_min;
Q = dat.cost.Q;
R = dat.cost.R;
P = dat.cost.P;
A = dat.model.A;
B = dat.model.B;
X = dat.track.x;
U = dat.track.u;
XINIT = dat.XINIT;
l = dat.l;
L = dat.L;

%% Plotting
figure;
plot(133+dxs(4,:),'r');
hold
plot(133+X(4,:));   
h1leg = legend('Reference', 'Generated trajectory');
h1title = title('Airspeed');
    
figure;
plot(dxs(7,:),'r');
hold
plot(X(7,:));
h1leg = legend('Reference', 'Generated trajectory');
h1title = title('Roll');



%% ----------------------------------------------
%% Solve w/ splitting methods - MPC for tracking
%% ----------------------------------------------
fcnList_MPC = { @templateAMA, @templateFAMA, @templatePrFAMA };

clear settings;
clear ITER_COUNT;

%% Iterate over the algorithms
for kkk = 1:length(fcnList_MPC)
    % Global settings for solver
    settings.maxItr = 1e3;
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
    
    %% Receding horizon control
    for ttt = 1:SimTime-2*N   
        if ttt == 1
            xinit = XINIT;
        end

        % Solve
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % FILL IN THE FUNCTIONS
        [sol, stats, stats_plot] =...
            fcnList_MPC{kkk}(dat,reshape(dxs(:,ttt:ttt+N),(N+1)*nx,1),...
            reshape(dus(:,ttt:ttt+N-1),N*nu,1),xinit,settings);
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        % Recover optimal closed loop trajectories
        fom.x{kkk} = reshape(sol.x(1:(N+1)*nx),nx,N+1);
        fama.x{kkk} = [fama.x{kkk} fom.x{kkk}(:,2)];
        fom.u{kkk} = reshape(sol.x((N+1)*nx+1:end),nu,N);
        fama.u{kkk} = [fama.u{kkk} fom.u{kkk}(:,1)];
        
        if (ttt==SimTime-2*N)
            fom.x{kkk} = reshape(sol.x(1:(N+1)*nx),nx,N+1);
            fama.x{kkk} = [fama.x{kkk} fom.x{kkk}(:,2:end)];
            fom.u{kkk} = reshape(sol.x((N+1)*nx+1:end),nu,N);
            fama.u{kkk} = [fama.u{kkk} fom.u{kkk}]; 
        end
        
        xinit = fom.x{kkk}(:,2)+1e-3*(-1+2*rand)*fom.x{kkk}(:,2); % perturbed initial state
        SOL{kkk,ttt} = sol;
        ITER_COUNT(kkk,ttt) = stats.numiter;
        STATS_PLOT{kkk,ttt,:} = stats_plot;
        ITER_COUNT(kkk,ttt)   = stats.numiter;
        DUAL_RESIDUAL(kkk,ttt) = sqrt(stats.rDual);
        ACTIVE_CONSTRAINTS(kkk,ttt) = sum(abs(sol.lam) > 1e-3);
    end
    %% Keep statistics for each algorithm
    ERROR(kkk) = norm([fama.x{kkk}] - X)  / norm(X);
    ITERS(kkk) = mean(ITER_COUNT(kkk,:));

    % forward simulate system with saturated (sub)optimal input sequence
    fama.u{kkk} = min(repmat(umax,1,195),max(fama.u{kkk},repmat(umin,1,195)));
    fama.x{kkk}(:,1) = XINIT;
    for ttt=1:195
        fama.x{kkk}(:,ttt+1) = A*fama.x{kkk}(:,ttt) + B*fama.u{kkk}(:,ttt);
    end

    fama.x{kkk}(4,:) = 133+fama.x{kkk}(4,:);
    fama.x{kkk}(7,:) = (180/pi)*fama.x{kkk}(7,:);
    fama.x{kkk}(8,:) = (180/pi)*fama.x{kkk}(8,:);
    Xorig(4,:) = 133+X(4,:);
    Xorig(7,:) = (180/pi)*X(7,:);
    Xorig(8,:) = (180/pi)*X(8,:);

    %% Plotting
    figure;
    plot(133+dxs(4,:),'r');
    hold
    plot(Xorig(4,:));   
    plot(fama.x{kkk}(4,:)); 
    h1leg = legend('Reference', 'IP solver', 'FAMA');
    h1title = title('Airspeed');

    figure;
    plot((180/pi)*dxs(7,:),'r');
    hold
    plot(Xorig(7,:));
    plot(fama.x{kkk}(7,:));
    h1leg = legend('Reference', 'IP solver', 'FAMA');
    h1title = title('Roll');
    pause;
end
fprintf('\n')
X_MPC = fom.x;
U_MPC = fom.u;
ERROR_MPC = ERROR;
ITERS_MPC = ITERS;

%% Plotting
figure(1);
hold on;
plot(1e3*X(7,:), 'k', 'LineWidth', 2);
plot(ITER_COUNT(1,:), '--g','LineWidth', 1.5);
plot(ITER_COUNT(2,:), '--b','LineWidth', 1.5);
plot(ITER_COUNT(3,:), '--r','LineWidth', 1.5);
hold off;
h1y = ylabel('No. of iterations'); h1x = xlabel('Time');
h1t = title('Iterations per problem');
h1leg = legend('Roll reference', 'AMA', 'FAMA', 'PrecFAMA');
set(h1t, 'FontSize', 14);