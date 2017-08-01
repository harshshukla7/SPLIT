% Generate custom c-code to solve ADMM
function sd = coder_CP(prob)

dat  = prob.coderData;
sd   = splitData;

%%%% First we find the constant values i.e. \rho \sigma \gamma 
L = []; l = []; 

prox = prob.prox;
nProx = length(prox);
kk = 1;
for i = 1:nProx
  l            = [l;prox(i).l];
  L            = [L;prox(i).L];
end

% define matrices specific for chambolle-pock
K_p_nu = [L; dat.A];
matrix_inverse = []; % CANNOT PRECOMPUTE INVERSE

% Bound on K operator. check CPI.m or CPII.m 
Mu = max(svds(K_p_nu)); 
tau = 1/Mu;
sigma = 1/(tau*Mu^2);
warning('Check the way gamma value is created says ''coder_CP''');
gamma = 0.1*min(eig(dat.Q));



a_name=sprintf('tau');
sd.add_var(a_name,tau,'type','real');
a_name=sprintf('sigma');
sd.add_var(a_name,sigma,'type','real');
a_name=sprintf('gamma');
sd.add_var(a_name,gamma,'type','real');

%%%%%%%% Next we find the dimensions


sd.define('nParam',  size(dat.pL,2), 'int');
sd.define('nPrimal', size(dat.A,2),  'int');
sd.define('nDual',   size(dat.L,1),  'int');
sd.define('nEqCon',  size(dat.A,1),  'int'); % this is nNu as well


%% This is to update p and nu. see step 1 of CPI.m. 
sd.add_function(coderFunc_times('custom_mult_K_p_nu', sparse(K_p_nu), ...
  'y_name', 'workDual', 'x_name', 'x', 'Astr', 'L'));

%%  Step 2 of CPI.m Evaluate prox functions y = prox(workDual)
sd.add_function(coderFunc_prox_CP(dat));

%% Step 3 of CPI.m Solve the KKT system

% Assumes that the KKT matrix is constant
%diagQ    = diag(dat.Q);
%diagItau = diag(eye(size(dat.Q, 1))*(1/tau));

diagItau = (eye(size(dat.Q, 1))*(1/tau));
K = (dat.Q + diagItau); %% 
%%%% calculate and transform the problem in suite sparse format.

LNZ=nnz(K)-nnz(triu(K));
ANZ=nnz(triu(K));
N=size(K,2);

%%%% add the varibles 

sd.define('nn_lp',size(K,2),'int');
sd.define('LNZ_ss', LNZ, 'int');
sd.define('ANZ_ss', ANZ, 'int');
sd.define('N_ss', N, 'int');

%  To compute K_p_nu'*x. See how x is computed in CPI.m
sd.add_function(...
  coderFunc_times('custom_mult_K_p_nu_trans', sparse(K_p_nu'), 'Astr', 'Ltrans', 'method', 'blas'));


sd.add_function(coderFunc_mldivide('custom_solve_kkt', K, 'Astr', 'KKT','method','ldl_ss'));

