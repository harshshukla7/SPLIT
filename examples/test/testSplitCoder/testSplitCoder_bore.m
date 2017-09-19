clear all

%%  Example 13.1 page 253

n = 2;
m = 1;
N = 3;

A = [1 1 ; 0 1];
B = [0;1];
C = [1,0];
D = 0;

x = splitvar(n,N);
u = splitvar(m,N-1);

%% Different constraints: 

%%%%% 1.  Box
-5 <= x(:) <= 5;
-0.5 <= u(:) <= 0.5;

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
dp = {'xo','tr_affine','tr_track'};

diff_parametrix = 'tr_track' ;

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
        
        objX = (x) .* (Q*(x));
        objU = u .* (R*u);
        obj = sum(objX(:)) + sum(objU(:));
        
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

% Generate the code into the splitData object sd
settings.adaptive = 'no';
settings.precond = 'no';
settings.E = eye(size(prob.coderData.L,1));
sd = coder_admm(prob,settings);

%sd_ama = coder_ama(prob);
%sd_CP = coder_CP(prob);
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
tic
sol    = AdPrADMM(prob);

%sol_f  = f_admm(prob);
%sol_ama = ama(prob);
%sol_fama = fama(prob);


%% For CPI and CPII

%sol_cpI = CPI(prob);
%sol_cpII = CPII(prob);

%%

sol_x  = sol.x;

sd.add_var('par_ex', par_ex);
sd.add_var('sol_x',  sol_x);

%%

% Write the generated files probData.c, probData.h and probData.dat
sd.write_to_file('probData');


% %% Compile command-line interface test file
% 
% % Compile command-line
% system('make clean; make')
% system('./admm')

%% Compile mex-interface

% Compile mex-interface

%mex -largeArrayDims admm_mex.c admm.c matrix_ops.c probData.c splitLoad.c splitTimer.c 
 mex -I/usr/local/include LDFLAGS='-L/usr/local/lib' -largeArrayDims admm_mex.c admm.c matrix_ops.c probData.c splitLoad.c splitTimer.c   -L/usr/local/lib ...
     LIBS='-lm -lldl -framework Accelerate '
% Test mex-interface

% opt.primalTol = 1e-4;
% opt.dualTol   = 1e-4;
% opt.MAXITR    = 1e3;
% opt.ITR_PER_CONV_TEST = 10;

% If you want to time just the solve-time, then set this larger than one,
% and the time returned will be the average per solve
% opt.number_of_solves = 100; 
% 
% 
% sol_mex = admm_mex(par_ex, opt)

%%
% fprintf('Error between c-code and m-code solutions : %e\n', norm(sol_mex.primal - sol.x));
% fprintf('Solve time matlab : %.2e us\n', solve_time_mat * 1e6);
% fprintf('Solve time split  : %.2e us\n', sol_mex.solve_time_ns / 1e3);
%%

% f = coderFunc('void myfunc(REAL x[%i])', 5)
% f.add_var('y', y); % Write the vector y to the file
% f.pl('for(int i=0; i<5; i++)')
% f.pl('  x[i] = y[i];')

