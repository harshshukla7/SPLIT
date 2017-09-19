clear all;  close all;
rand('state', 0);
randn('state', 0);

fcnList_MPC = { @PrFAMA, @AdPrADMM, @AdPrFADMM, @FPDA};

countVars = 1;
for iii = 30%2:2:30 % nx+nu
    countHorzs = 1;
    for ttt = 26%2:2:30 % N
        nx = iii;
        nu = 2*nx;
        sys{iii,ttt} = drss(nx,nx,nu);
        A{iii,ttt} = randn(2*(nx+nu),(nx+nu));
        b{iii,ttt} = 0.1*iii*abs(randn(2*(nx+nu),1));
        alpha{iii,ttt} = [abs(randn);abs(randn)];
        z0{iii,ttt} = randn(nx,1);
%         randQ = randn(nx,nx);
%         randR = randn(nu,nu);
        Q{iii,ttt} = eye(nx);   R{iii,ttt} = eye(nu);


        %% Solve CVX
        cvx_begin
        cvx_solver mosek
        variables X(nx,ttt+1) U(nu,ttt)
        obj = 0;
        for j = 1:ttt
           obj = obj + 0.5*X(:,j)'*Q{iii,ttt}*X(:,j)+0.5*U(:,j)'*R{iii,ttt}*U(:,j)+alpha{iii,ttt}(1)*norm(U(:,j),1)+alpha{iii,ttt}(2)*norm(X(:,j),2);
        end
        obj = obj + 0.5*X(:,ttt+1)'*Q{iii,ttt}*X(:,ttt+1)+alpha{iii,ttt}(2)*norm(X(:,ttt+1),2);
        minimize(obj)
        subject to
        X(:,1) == z0{iii,ttt};
        for j = 1:ttt
            X(:,j+1) == sys{iii,ttt}.A*X(:,j) + sys{iii,ttt}.B*U(:,j);
            U(:,j) <= iii*0.1;
            U(:,j) >= -iii*0.1;
            norm(X(:,j+1),2) <= 2;
        end
        cvx_end
            
        %% Solve with FoMs
        % Set the problem in split format
        for kkk = 1:length(fcnList_MPC)         
            % Feasibility check
            if ~(strcmp(cvx_status,'Solved') || strcmp(cvx_status,'Inaccurate/Solved'))
                ITER_COUNT(:,countVars,countHorzs) = 0;
                break;
            end
            
            splitProb.clearProblem;
            x = splitvar(nx,ttt+1);
            u = splitvar(nu,ttt);  
            obj = 0;
            for j = 1:ttt
                obj = obj + 0.5*x(:,j)'*Q{iii,ttt}*x(:,j)+0.5*u(:,j)'*R{iii,ttt}*u(:,j)+alpha{iii,ttt}(1)*norm(u(:,j),1)+alpha{iii,ttt}(2)*norm(x(:,j),2);
            end
            obj = obj + 0.5*x(:,ttt+1)'*Q{iii,ttt}*x(:,ttt+1)+alpha{iii,ttt}(2)*norm(x(:,ttt+1),2);
            x(:,1) == z0{iii,ttt};
            for j = 1:ttt
                x(:,j+1) == sys{iii,ttt}.A*x(:,j) + sys{iii,ttt}.B*u(:,j);
                u(:,j) <= iii*0.1;
                u(:,j) >= -iii*0.1;
                norm(x(:,j+1),2) <= 2;
            end
            minimize(obj);
            prob = splitProb.genProblem;
            settings.maxItr  = 2e3;
            settings.dualTol = 1e-2;
            settings.primalTol = 1e-2;
            
            %% settings
            if (kkk==1) % FAMA
                settings.restart = 'yes';
            %     settings = rmfield(settings, 'restart'); 
            elseif (kkk==2) % ADMM
                settings.relaxation  = 1.8;
                settings.rho = 1;
            %         settings = rmfield(settings, 'adapt');
            elseif (kkk==3) % FADMM
                settings.rho = 1;
            elseif (kkk==4) % FPDA
                settings = rmfield(settings, 'rho');
            end
            
            % Solve
            [sol, stats, stats_plot] = fcnList_MPC{kkk}(prob, settings);
            
            % Stats
            fom.x{kkk,countVars,countHorzs} = reshape(sol.x(1:(ttt+1)*nx,1),nx,ttt+1); 
            fom.u{kkk,countVars,countHorzs} = reshape(sol.x(nx*(ttt+1)+1:nx*(ttt+1)+ttt*nu,1),nu,ttt);
            ITER_COUNT(kkk,countVars,countHorzs) = stats.numiter;
            STATS_PLOT{kkk,countVars,countHorzs} = stats_plot;
            ERROR(kkk,countVars,countHorzs) = (norm([fom.u{kkk,countVars,countHorzs}] - U) + norm([fom.x{kkk,countVars,countHorzs}] - X))/(norm(U)+norm(X));
%             ITERS(kkk,iii,ttt) = mean(ITER_COUNT(kkk,:,:));
        end
        countHorzs = countHorzs+1; 
    end
    countVars = countVars+1; 
end

%% plotting
Horzs = [2:2:30];
Vars = [2:2:30];
[X,Y] = meshgrid(Horzs, Vars);
Titles = {'FAMA', 'ADMM', 'FADMM', 'FPDA'};
Z = {};

for iii = 1:4
    figure;
    Z{iii} = squeeze(ITER_COUNT(iii,:,:))';
    M(iii) = median(median(Z{iii}));
    surf(X,Y,Z{iii})
    xlabel('Horizon')
    ylabel('Variables')
    colormap default
    title(Titles{iii});
    caxis([9 2e3])
end