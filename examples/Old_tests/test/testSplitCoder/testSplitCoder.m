clear all

dp = {'x0','tr_affine','tr_track'};

diff_parametrix = 'x0' ;

algorithms = {'ADMM','FADMM','FAMA','PDA','FPDA'};

swit_algo = 'ADMM';

compile_using = {'system','mex','none'};

compile_option = 'none';

n = 5;
m = 3;
N = 10;

[A,B,C,D] = ssdata(drss(n,n,m));

x = splitvar(n,N);
u = splitvar(m,N-1);

%% Different constraints:

%%%%% 1.  Box
-10 <= x(:) <= 10;
-1 <= u(:) <= 1;

%%%%% 2. l1-ball
%F = 9e-3*randn(n*N,n);
%f = 1e-6*randn(n,1);
%d = 100;
%norm(F'*x(:)+f,1) <= d;

%%%%% 2. l2-ball
%F = 50*randn(n*N,n);
%f = 1e-7*randn(n,1);
%d = 1;
%norm(F'*x(:)+f,2) <= d;

%%%%% 3. linf ball
%F = 15*randn(n*N,n);
%f = 1e-5*randn(n,1);
%d = 1;
%norm(F'*x(:)+f,inf) <= d;

%%%%% 4. Lorentz Cone

%F = 50*randn(n*N,n);
%f = 1e-7*randn(n,1);
%d = 1;
%q = 1e-3*randn(n*N,1);

%norm(F'*x(:)+f,2) <= q'*x(:) + d;


%%%%% 5. Polytope
%K = 1e-1*randn(3*(n*N+m*(N-1)),n*N+m*(N-1));
%k = 20*abs(randn(3*(n*N+m*(N-1)),1));
%K*[x(:); u(:)] <= k;

%%


switch diff_parametrix
    case 'x0'
        
        
        x0 = parameter(n,1);
        
        %r = parameter(n,N);
        x(:,1) == x0;
        x(:,2:end) == A*x(:,1:end-1) + B*u;
        
        Q = eye(n);
        R = eye(m);
        % objX = x .* (Q*x);
        % objU = u .* (R*u);
        % obj = sum(objX(:)) + sum(objU(:));
        
%         objX = x.*(Q*(x));
%         objU = u.*((R*u));
%         obj = sum(objX(:)) + sum(objU(:));
obj = 0;
for i=1:N-1
    obj = obj + x(:,i)'*Q*x(:,i) + u(:,i)'*R*u(:,i);
end
obj = obj + x(:,N)'*Q*x(:,N);

        
        minimize(obj);
        
        
    case 'tr_affine'
        
        
        r = parameter(n,N);
        
        Q = eye(n);
        R = eye(m);
        
        objX = (x) .* (Q*(x));
        objU = u .* (R*u);
        objA = r .* x;
        obj =  + sum(objA(:));
        
        minimize(obj);
        
    case 'tr_track'
        
        x0 = parameter(n,(N+1)); % first n for initia state and the last for reference signal
        
        %r = parameter(n,N);
        x(:,1) == x0(:,1);
        x(:,2:end) == A*x(:,1:end-1) + B*u;
        
        Q = eye(n);
        R = eye(m);
        % objX = x .* (Q*x);
        % objU = u .* (R*u);
        % obj = sum(objX(:)) + sum(objU(:));
        
        objX = (x) .* (Q*(x));
        objU = u .* (R*u);
        objPar = x0(:,2:(N+1)).* ((2*Q)*x);
        obj = sum(objX(:)) + sum(objU(:)) - sum(objPar(:));
        
        minimize(obj);
        
end


prob = splitProb.genProblem;


%%

% Add a particular problem that we're going to solve

switch diff_parametrix
    case 'x0'
        par_ex = 10*ones(n,1);
        x0.set(par_ex);
        
    case 'tr_affine'
        par_ex = 10*ones(n,N);
        r.set(par_ex);
        
    case 'tr_track'
        par_ex = zeros(n,N+1);
        par_ex(:,1) = 1*ones(n,1);
        par_ex(:,2:end) = 5*ones(n,N);
        x0.set(par_ex);
end

%%

switch swit_algo
    
    case 'ADMM'
        settings.adaptive = 'no';
        settings.precond = 'no';
        settings.E = eye(size(prob.coderData.L,1));
        sd = coder_admm(prob,settings);
        sol    = AdPrADMM(prob);
        sol_x  = sol.x;
        
        sd.add_var('par_ex', par_ex);
        sd.add_var('sol_x',  sol_x);
        
        % Write the generated files probData.c, probData.h and probData.dat
        sd.write_to_file('probData');
        
    case 'FADMM'
        settings.adaptive = 'no';
        settings.precond = 'no';
        
        settings.E = eye(size(prob.coderData.L,1));
        sd = coder_admm(prob,settings);
        sol    = AdPrFADMM(prob);
        sol_x  = sol.x;
        
        sd.add_var('par_ex', par_ex);
        sd.add_var('sol_x',  sol_x);
        
        % Write the generated files probData.c, probData.h and probData.dat
        sd.write_to_file('probData');
        
    case 'FAMA'
        tp_settings.precond = 'no';
        tp_settings.Mat_Vec = 'ss';
        tp_settings.Lin_Solve = 'ldl_lp';
        sd_ama = coder_ama(prob,tp_settings);
        set_ama.rho = 1;
        sol_ama = PrFAMA(prob, set_ama);
        sol  = sol_ama;
        
        sol_ama_x = sol_ama.x;
        sd_ama.add_var('par_ex', par_ex);
        sd_ama.add_var('sol_x',  sol_ama_x);
        
        sd_ama.write_to_file('probData');
        
    case 'PDA'
        clear tp_settings;
        tp_settings.precond = 'no';
        tp_settings.Mat_Vec = 'ss';
        tp_settings.Lin_Solve = 'ldl_lp';
        set_PDA = [];
        sol_PDA= AdPrPDA(prob, set_PDA);
        
        sol = sol_PDA;
         clear sd_PDA
         sd_PDA = coder_PDA(prob,tp_settings);
         sol_PDA_x = sol_PDA.x;
         sd_PDA.add_var('par_ex', par_ex);
         sd_PDA.add_var('sol_x',  sol_PDA_x);
%         
         sd_PDA.write_to_file('probData');
   
    case 'FPDA'
        clear tp_settings;
        tp_settings.precond = 'no';
        tp_settings.Mat_Vec = 'ss';
        tp_settings.Lin_Solve = 'ldl_lp';
        set_FPDA = [];
        sol_FPDA= FPDA(prob, set_FPDA);
        sol = sol_FPDA;
         clear sd_FPDA
         sd_FPDA = coder_FPDA(prob,tp_settings);
         sol_FPDA_x = sol_FPDA.x;
         sd_FPDA.add_var('par_ex', par_ex);
         sd_FPDA.add_var('sol_x',  sol_FPDA_x);
%         
         sd_FPDA.write_to_file('probData');
   
end


switch compile_option
    
    case 'none'
        
    case 'system'
        %% Compile command-line interface test file
        
        % Compile command-line
        system('make clean; make')
        system('./fama')
        
        
    case 'mex'
        %% Compile mex-interface
        switch swit_algo
            
            case 'ADMM'
                
%       mex -Iinclude/ilcplex cplexint.c lib/x86-64_sles10_4.1/static_pic/libcplex.a         

%                 mex  -v -I/usr/local/include LDFLAGS='-L/usr/local/lib'  -largeArrayDims admm_mex.c admm.c matrix_ops.c probData.c splitLoad.c splitTimer.c   -L/usr/local/lib ...
%      LIBS='-lm -lldl -framework Accelerate' -lmwlapack -lmwblas

        mex  -v -I/usr/local/include LDFLAGS='-L/usr/local/lib'  -largeArrayDims admm_mex.c admm.c matrix_ops.c probData.c splitLoad.c splitTimer.c   -L/usr/local/lib ...
                -lm -lldl -lmwlapack -lmwblas
                
               % Test mex-interface
                
                opt.primalTol = 1e-4;
                opt.dualTol   = 1e-4;
                opt.MAXITR    = 1e3;
                opt.ITR_PER_CONV_TEST = 10;
                
                %If you want to time just the solve-time, then set this larger than one,
                %and the time returned will be the average per solve
                opt.number_of_solves = 100;
                
                
                sol_mex = admm_mex(par_ex, opt);
                
            case 'FADMM'
               
                mex admm_mex -largeArrayDims admm_mex.c fadmm.c matrix_ops.c probData.c splitLoad.c splitTimer.c
                
               % Test mex-interface
                
                opt.primalTol = 1e-4;
                opt.dualTol   = 1e-4;
                opt.MAXITR    = 1e3;
                opt.ITR_PER_CONV_TEST = 10;
                
                %If you want to time just the solve-time, then set this larger than one,
                %and the time returned will be the average per solve
                opt.number_of_solves = 100;
                
                
                sol_mex = admm_mex(par_ex, opt);
            case 'FAMA'
                
                mex  -v  -largeArrayDims ama_mex.c fama.c matrix_ops.c probData.c splitLoad.c splitTimer.c   -L/usr/local/lib ...
     LIBS='-lm -lldl -framework Accelerate' -lmwlapack -lmwblas -I/usr/local/include LDFLAGS='-L/usr/local/lib'

                
                Test mex-interface
                
                opt.primalTol = 1e-4;
                opt.dualTol   = 1e-4;
                opt.MAXITR    = 1e3;
                opt.ITR_PER_CONV_TEST = 10;
                
                %If you want to time just the solve-time, then set this larger than one,
                %and the time returned will be the average per solve
                opt.number_of_solves = 100;
                
                
                sol_mex = ama_mex(par_ex, opt);
                
            case 'PDA'
                
                
                opt.primalTol = 1e-3;
                opt.dualTol   = 1e-3;
                opt.MAXITR    = 500;
                opt.ITR_PER_CONV_TEST = 5;
                
                %If you want to time just the solve-time, then set this larger than one,
                %and the time returned will be the average per solve
                opt.number_of_solves = 1;
                
                mex -output PDA_mex  -v  -largeArrayDims  PDA.c matrix_ops.c probData.c  splitLoad.c splitTimer.c ama_mex.c ldl.c

                sol_mex = PDA_mex(par_ex, opt);
                

            case 'FPDA'
                opt.primalTol = 1e-3;
                opt.dualTol   = 1e-3;
                opt.MAXITR    = 500;
                opt.ITR_PER_CONV_TEST = 5;
                
                %If you want to time just the solve-time, then set this larger than one,
                %and the time returned will be the average per solve
                opt.number_of_solves = 1;
                
                mex -output FPDA_mex  -v  -largeArrayDims  FPDA.c matrix_ops.c probData.c  splitLoad.c splitTimer.c ama_mex.c ldl.c

                sol_mex = FPDA_mex(par_ex, opt);
                
        end
        %%
            fprintf('Error between c-code and m-code solutions : %e\n', norm(sol_mex.primal - sol.x));
%            fprintf('Solve time matlab : %.2e us\n', solve_time_mat * 1e6);
            fprintf('Solve time split  : %.2e ms\n', sol_mex.solve_time_ns / 1e6);
        %%
        
        % f = coderFunc('void myfunc(REAL x[%i])', 5)
        % f.add_var('y', y); % Write the vector y to the file
        % f.pl('for(int i=0; i<5; i++)')
        % f.pl('  x[i] = y[i];')
        
        
        
end
