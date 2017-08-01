%% to compile c code following changes need to be done
%%%%%%%%%%%%%%%%%%%%%% MAKE FILE %%%%%%%%%%%%%%%%%%%%%%
% 1. Line 37: change the .c according to FPDA,FADMM etc
% 2. Line 50: Name
% 3. Line 61: Echo
%%%%%%%%%%%%%%%%%%%%%% MAIN %%%%%%%%%%%%%%%%%%%%%%
% 1. Line 83: kkk value
% 2. Line 159: ./FPDA etc
% 3. Line 121: change the method for the linear solve step
% 4. Line 20: change new dimension value = 1 for new simulations (for new)
% 5. Line 29: change state dimensions (for a new simulations)
% 6. Line 32: change horizion length  (for a new simulations)
%%%%%%%%%%%%%%%%%%%%%% test_coder.c %%%%%%%%%%%%%%%%%%%%%%
% 1. line 3: .h file inclusion
% 2. line 70: change the name of the text file written
% 3. line 16: change the number of simulation iterations

%%%%%% change the line 28

new_dimension = 1;

if new_dimension
    clear all;  close all;
    rand('state', 0);
    randn('state', 0);
    
    cvx_do = 1;
    fcnList_MPC = { @PrFAMA, @AdPrADMM, @AdPrFADMM, @FPDA};
    
    countVars = 1;
    iii =  2;  %: 5 : 20 % nx
    countHorzs = 1;
    ttt =  7; %: 5 : 20 % N 30
    nx = iii;
    %nu = 2*nx;
    nu = nx;
    nx = 2*nu;
    bd = nx;
    bd_u = nu;
    
    sys{iii,ttt} = drss(nx,nu,nu);
    A{iii,ttt} = randn(2*(nx+nu),(nx+nu));
    b{iii,ttt} = 0.1*iii*abs(randn(2*(nx+nu),1));
    alpha{iii,ttt} =  [abs(randn);abs(randn)];
    %z0{iii,ttt} = 200*randn(nx,1);
    z0{iii,ttt} = randn(nx,1);
    %         randQ = randn(nx,nx);
    %         randR = randn(nu,nu);
    Q{iii,ttt} = eye(nx);
    R{iii,ttt} = eye(nu);
    %FF = 1e-3*createbandOPT(nx,bd);
    %FF_u = 1e-3*createbandOPT(nu,bd_u);
    FF = 1e-5*eye(nx,nx);
    FF_u = 1e-5*eye(nu,nu);
    FF_x_hs = 0.01;
    FF_u_hs = 0.01;
    hs_scale = 0.01;
    
    %% Solve CVX
    
    if(cvx_do == 1)
        cvx_begin
        cvx_solver sdpt3
        variables X(nx,ttt+1) U(nu,ttt)
        obj = 0;
        for j = 1:ttt
            obj = obj + 0.5*X(:,j)'*Q{iii,ttt}*X(:,j)+0.5*U(:,j)'*R{iii,ttt}*U(:,j); %+alpha{iii,ttt}(1)*norm(FF_u*U(:,j),1)+alpha{iii,ttt}(2)*norm(X(:,j),2);
        end
        obj = obj + 0.5*X(:,ttt+1)'*Q{iii,ttt}*X(:,ttt+1); %+alpha{iii,ttt}(2)*norm(X(:,ttt+1),2);
        minimize(obj)
        subject to
        X(:,1) == z0{iii,ttt};
        for j = 1:ttt
            X(:,j+1) == sys{iii,ttt}.A*X(:,j) + sys{iii,ttt}.B*U(:,j);
            U(:,j) <= iii*hs_scale;
            U(:,j) >= -iii*hs_scale;
            % norm(FF*X(:,j+1),2) <= FF_x_hs;
            % norm(FF_u*U(:,j),2) <= FF_u_hs;
        end
        cvx_end
    end
    
end
%
%% Solve with FoMs
% Set the problem in split format
kkk = 1; %:length(fcnList_MPC) { @PrFAMA, @AdPrADMM, @AdPrFADMM, @FPDA};
% Feasibility check
%             if ~(strcmp(cvx_status,'Solved') || strcmp(cvx_status,'Inaccurate/Solved'))
%                 ITER_COUNT(:,countVars,countHorzs) = 0;
%                 break;
%             end

splitProb.clearProblem;
x = splitvar(nx,ttt+1);
u = splitvar(nu,ttt);
x_par = parameter(nx,1);
obj = 0;
for j = 1:ttt
    obj = obj + 0.5*x(:,j)'*Q{iii,ttt}*x(:,j)+0.5*u(:,j)'*R{iii,ttt}*u(:,j); %   + alpha{iii,ttt}(1)*norm(u(:,j),1)+alpha{iii,ttt}(2)*norm(x(:,j),2);
end
obj = obj + 0.5*x(:,ttt+1)'*Q{iii,ttt}*x(:,ttt+1);  %+alpha{iii,ttt}(2)*norm(x(:,ttt+1),2);
x(:,1) == x_par;
for j = 1:ttt
    x(:,j+1) == sys{iii,ttt}.A*x(:,j) + sys{iii,ttt}.B*u(:,j);
    u(:,j) <= iii*hs_scale;
    u(:,j) >= -iii*hs_scale;
    %                 norm(FF*x(:,j+1),2) <= FF_x_hs;
    %norm(FF_u*u(:,j),2) <= FF_u_hs;
    %                norm(x(:,j+1),2) <= 2;
    %               norm(u(:,j),2) <= 2;
end
minimize(obj);
prob = splitProb.genProblem;
settings.maxItr  = 2e3;
settings.dualTol = 1e-2;
settings.primalTol = 1e-2;
x_par.set(z0{iii,ttt});
%%% settings
if (kkk==1) % FAMA
    %settings.restart = 'yes';
    %     settings = rmfield(settings, 'restart');
    settings.precond = 'no';
    settings.Mat_Vec = 'ss'; % FPGA_matvec ss
    settings.Lin_Solve = 'invert_FPGA'; % ldl_ss, invert, ldl, ldl_lp, invert_FPGA
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
sd.add_var('par_ex',z0{iii,ttt});
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
fom.x{kkk,countVars,countHorzs} = reshape(sol.x(1:(ttt+1)*nx,1),nx,ttt+1);
fom.u{kkk,countVars,countHorzs} = reshape(sol.x(nx*(ttt+1)+1:nx*(ttt+1)+ttt*nu,1),nu,ttt);
ITER_COUNT(kkk,countVars,countHorzs) = stats.numiter;
STATS_PLOT{kkk,countVars,countHorzs} = stats_plot;
ERROR(kkk,countVars,countHorzs) = (norm([fom.u{kkk,countVars,countHorzs}] - U) + norm([fom.x{kkk,countVars,countHorzs}] - X))/(norm(U)+norm(X));
%             ITERS(kkk,iii,ttt) = mean(ITER_COUNT(kkk,:,:));

countHorzs = countHorzs+1;

countVars = countVars+1;

