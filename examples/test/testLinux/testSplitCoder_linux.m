clear all

n = 3;
m = 3;
N = 10;

[A,B,C,D] = ssdata(drss(n,n,m));

x = splitvar(n,N);
u = splitvar(m,N-1);
x0 = parameter(n,1);

x(:,1) == x0;

%% Different constraints: 

%%%%% 1.  Box
%-10 <= x(:) <= 10;
%-1 <= u(:) <= 1;

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
K = 1e-1*randn(3*(n*N+m*(N-1)),n*N+m*(N-1));
k = 20*abs(randn(3*(n*N+m*(N-1)),1));
K*[x(:); u(:)] <= k;

%%
x(:,2:end) == A*x(:,1:end-1) + B*u;

Q = eye(n);
R = eye(m);
objX = x .* (Q*x);
objU = u .* (R*u);
obj = sum(objX(:)) + sum(objU(:));

minimize(obj);

prob = splitProb.genProblem;




%%

settings.precond = ('no');
settings.adaptive = ('no');
% Generate the code into the splitData object sd
sd = coder_admm(prob,settings);

%sd_ama = coder_ama(prob);
%sd_CP = coder_CP(prob);
%%

% Add a particular problem that we're going to solve
par_ex = 1*ones(n,1);
x0.set(par_ex);
sol    = admm(prob);
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


%%

% f = coderFunc('void myfunc(REAL x[%i])', 5)
% f.add_var('y', y); % Write the vector y to the file
% f.pl('for(int i=0; i<5; i++)')
% f.pl('  x[i] = y[i];')