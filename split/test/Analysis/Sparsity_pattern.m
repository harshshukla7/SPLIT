clear all

n = 5;
m = 5;
N = 15;

[A,B,C,D] = ssdata(drss(n,n,m));

x = splitvar(n,N);
u = splitvar(m,N-1);
x0 = parameter(n,1);

x(:,1) == x0;

%% Different constraints: 

 constra = 'box'; % 'box', 'l1', 'l2', 'linf', 'lorentz', 'polytope', 'polytope_sparse'

 switch constra
     
%%%%% 1.  Box
    case 'box'
-10 <= x(:) <= 10;
-1 <= u(:) <= 1;

 
%%%%% 2. l1-ball
     case 'l1'
F = 9e-3*randn(n*N,n);
f = 1e-6*randn(n,1);
d = 100;
norm(F'*x(:)+f,1) <= d;

%%%%% 2. l2-ball
     case 'l2'
F = 50*randn(n*N,n);
f = 1e-7*randn(n,1);
d = 1;
norm(F'*x(:)+f,2) <= d;

%%%%% 3. linf ball
     case 'linf'
F = 15*randn(n*N,n);
f = 1e-5*randn(n,1);
d = 1;
norm(F'*x(:)+f,inf) <= d;

%%%%% 4. Lorentz Cone
     case 'lorentz'
F = 50*randn(n*N,n);
f = 1e-7*randn(n,1);
d = 1;
q = 1e-3*randn(n*N,1);

norm(F'*x(:)+f,2) <= q'*x(:) + d;


%%%%% 5. Polytope
     case 'polytope'
K = 1e-1*randn(3*(n*N+m*(N-1)),n*N+m*(N-1));
k = 20*abs(randn(3*(n*N+m*(N-1)),1));
K*[x(:); u(:)] <= k;

     case 'polytope_sparse'

K = 1e-1*eye(3*(n*N+m*(N-1)),n*N+m*(N-1));
k = 20*abs(randn(3*(n*N+m*(N-1)),1));
K*[x(:); u(:)] <= k;

 end
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

% Generate the code into the splitData object sd
sd = coder_admm(prob);

%sd_ama = coder_ama(prob);
%sd_CP = coder_CP(prob);
%%

% Add a particular problem that we're going to solve
par_ex = 10*ones(n,1);
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


%% Matrix - vect multiplication

dat  = prob.coderData;

figure()


subplot(1,2,1)
spy(dat.L');
title('dat\_L');
subplot(1,2,2)
spy(dat.L);
title('dat\_L');

%% KKT Analysis

rho = 1.0;
Q   = dat.Q;
rho = rho;
L   = dat.L;
A   = dat.A;
K = [Q+rho*L'*L A'; A zeros(size(A,1))];
[Kl, Kll, Klll] = ldl(K);
K1 = [Q A'; A zeros(size(A,1))];
[K1l, K1ll, K1lll] = ldl(K1);

figure()
subplot(1,2,1)
spy(K)
title('K\_admm');
subplot(1,2,2)
spy(Kl)
title('its L decomposition after permutation')

figure()
subplot(1,2,1)
spy(K1)
title('K\_ama');
subplot(1,2,2)
spy(K1l)
title('its L decomposition after permutation')
%% Reordering of the matrix : 

pk = symrcm(K);

figure()
subplot(1,2,1)
spy(K)
title('K matrix')
subplot(1,2,2)
spy(K(pk,pk))
title('K after permutation')
% f = coderFunc('void myfunc(REAL x[%i])', 5)
% f.add_var('y', y); % Write the vector y to the file
% f.pl('for(int i=0; i<5; i++)')
% f.pl('  x[i] = y[i];')