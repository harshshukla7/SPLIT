function [sol, stats_out, stats_plot, data_to_export] = ama_old(prob, settings)
%
% Quick implementation of ADMM to test the export from the split user
% interface
%
% Solves the problem
%   min 0.5 x' * Q * x + f' * x + sum_i Fi(Li * x + li)
%   s.t.  A * x = b
%
% prob.dat   : Contains Q, c, A and b matrices
% prob.prox  : Array of structures
%    prox.L, l Data matrices
%    prox.func - function handle to proximal function
%
% settings.rho = stepsize
%
% If the problem is parametric, then the current values of the parameters
% are used when solving the problem. Use set to set the parameter values
%

if nargin < 2
  settings = [];    
end

% The problem is parametric
if isfield(prob, 'paramInd') && length(prob.paramInd) > 0
  prob = setProbParameters(prob);
end

%%% NOTE: setting relaxation to something other than 1 seems to make SDPs fail

if ~isfield(settings, 'terminationCheckFreq'), settings.terminationCheckFreq = 10; end
if ~isfield(settings, 'rho'),                  settings.rho = 1; end
if ~isfield(settings, 'relaxation'),           settings.relaxation = 1; end
if ~isfield(settings, 'dualTol'),              settings.dualTol = 1e-3; end
if ~isfield(settings, 'primalTol'),            settings.primalTol = 1e-3; end
if ~isfield(settings, 'maxItr'),               settings.maxItr = 1e3; end

% Convergence statisctics requested
RECORD_STATS = false;
if nargout > 1 
  RECORD_STATS = true;
end

dat   = prob.dat;
prox  = prob.prox;
nProx = length(prox);

% Validate inputs
validate_inputs(prob, settings);

% Pull out parameters
% rho   = settings.rho;
% alpha = settings.relaxation;

% -----------------------
% Pre-processing
% -----------------------

n = size(dat.Q,1);
m = size(dat.A,1);

% Setup the KKT matrix and pre-solve
Q = dat.Q;
% sum_l = zeros(n,1);
sum_L = [];
L = []; l = []; 
nProxVars = zeros(nProx,1);
kk = 1;

for i = 1:nProx
  %Q            = Q + rho*prox(i).L'*prox(i).L;
  %sum_l        = sum_l + rho*prox(i).L'*prox(i).l;
  sum_L        = [sum_L prox(i).L'];
  l            = [l;prox(i).l];
  L            = [L;prox(i).L];
  proxInd{i}   = kk + [0:size(prox(i).L,1)-1];
  kk           = kk + size(prox(i).L,1);
  nProxVars(i) = size(prox(i).L,1);
end
proxWeight = [prob.prox.weight]; % Weights for the prox functions

% Compute the best stepsize:
eigenvalue_Q = eig(Q);
eigenvalue_LL = eig(L'*L);
% rho = = eigenvalue_Q(1)/eigenvalue_LL(end);

% Test: tune the stepsize rho:
rho = eigenvalue_Q(1)/eigenvalue_LL(end);
settings.rho = rho;

% Solution to KKT system:
%   [x;dual] = iKKT * (rhsKKTconst + rhsKKTvar * q)
% where q = lambda - y
KKT         = [Q dat.A'; dat.A zeros(m)];
inv_KKT = inv(KKT);
rhsKKTconst = [-dat.f' ;dat.b];
rhsKKTvar   = [sum_L;sparse(length(dat.b),size(sum_L,2))];

% Form the augmented lagrangian
% min 0.5*x'*Q*x + f'*x + sum_i Fi(y_i) + rho/2*sum_i ||Li*x + li - y_i + lam_i||
%      A*x = b

% Matrices to store the solution
x       = zeros(n,1);
y       = zeros(sum(nProxVars),1);
lam     = zeros(sum(nProxVars),1);
q       = zeros(size(L,1),1);

% precompute some common expression
scaled_L     = -rho*sum_L;
const_vector = inv_KKT * rhsKKTconst;
matrix       = inv_KKT * rhsKKTvar;

iii = 0;
for jjj = 1:settings.maxItr
% for jjj = 1:1
  iii = iii + 1;
  
  prev_y = y;

  % Step 1 : solve linear system
  % modification by Alex: I am not pretty sure, whether it is the best way
  % to do it - ask the guys !!!!
  % I have expanded Colin's formulation to prepare data for an export to
  % a header file
  t = const_vector + matrix * lam;
  x = t(1:n);
  
  % Compute over-relaxation
  %Lx_hat = alpha*L*x + (1 - alpha)*y;
    
  % Step 2 : compute prox operators
  q = L*x + l - 1/rho * lam;
  
  for i = 1:nProx
    ind = proxInd{i};
    % y(ind) = prox(i).func(q(ind)); 
    y(ind) = prox(i).func(q(ind), 1/rho*proxWeight(i));
  end
  
  % Update the lagrange multiplier
  prev_lam = lam;
  lam = rho*(y - q); 
    
  % Convergence checks
  s = SP*sum_L*(prev_lam - lam);
  r = invSD*(L*x + l - y);
  
  rPrimal = norm(r);
  rDual   = norm(s);

  % we need record stats inside the loop as well for plots
  if RECORD_STATS
     stats_plot(iii).rDual   = rDual;
     stats_plot(iii).rPrimal = rPrimal;
     stats_plot(iii).x       = x;
     stats_plot(iii).y       = y;
     stats_plot(iii).lam     = lam;
     stats_plot(iii).s       = s;
     stats_plot(iii).r       = r;
     stats_plot(iii).rho     = rho;
  end
  
  if rPrimal < settings.primalTol && rDual < settings.dualTol
    cprintf([0.1 0.5 0.1], '\n\n>>>>> Stopping on optimality after %i steps <<<<<\n\n', iii);
    break
  end
 
end

% we need record stats inside the loop as well for plots
if RECORD_STATS
   stats_out.rDual   = rDual;
   stats_out.rPrimal = rPrimal;
   stats_out.x       = x;
   stats_out.y       = y;
   stats_out.lam     = lam;
   stats_out.s       = s;
   stats_out.r       = r;
   stats_out.rho     = rho;
   stats_out.numiter = jjj;
end

sol.problem = 0;
if rPrimal > settings.primalTol || rDual > settings.dualTol
  cprintf('err', 'Failed to reach optimality after %i iterations\n', settings.maxItr);
  sol.problem = 1;
end

sol.lam = lam;
sol.x = x;
sol.y = y;

% prepare data for export to C

% common data required to solve the problem, e.g. stepsize, number of
% states, inputs ...
data_to_export.num_of_opt_var = n;
data_to_export.num_of_states  = m;
data_to_export.nProx          = nProx;
data_to_export.settings       = settings;

% sizes of vectors
data_to_export.b_LEN          = length(dat.b);
data_to_export.f_LEN          = length(dat.f);
data_to_export.c_LEN          = length(dat.c);
data_to_export.l_LEN          = length(l);
data_to_export.x_LEN          = length(x);
data_to_export.y_LEN          = length(y);
data_to_export.prev_y_LEN     = length(y);
data_to_export.lam_LEN        = length(lam);
data_to_export.q_LEN          = length(q);
data_to_export.t_LEN          = length(t);
data_to_export.const_vec_LEN  = length(const_vector);

% sizes of matrices
data_to_export.A_ROWS         = size(dat.A, 1);
data_to_export.A_COLS         = size(dat.A, 2);
data_to_export.Q_ROWS         = size(dat.Q, 1);
data_to_export.Q_COLS         = size(dat.Q, 2);
data_to_export.L_ROWS         = size(L, 1);
data_to_export.L_COLS         = size(L, 2);
data_to_export.scaledL_ROWS   = size(scaled_L, 1);
data_to_export.scaledL_COLS   = size(scaled_L, 2);
data_to_export.matrix_ROWS    = size(matrix, 1);
data_to_export.matrix_COLS    = size(matrix, 2);

% vectors
data_to_export.x              = zeros(n,1);
data_to_export.y              = zeros(sum(nProxVars),1);
data_to_export.prev_y         = zeros(sum(nProxVars),1);
data_to_export.lam            = zeros(sum(nProxVars),1);
data_to_export.q              = zeros(size(L,1),1);
data_to_export.b              = full(dat.b);
data_to_export.f              = full(dat.f);
data_to_export.c              = full(dat.c);
data_to_export.l              = full(l);
data_to_export.t              = zeros(1, data_to_export.t_LEN);
data_to_export.const_vec      = full(const_vector);

new_prox_beg   = [];
new_prox_end   = [];
len_of_vectors = [];
temp_sum = 0;
for i = 1:nProx
    new_prox_beg    = [ new_prox_beg proxInd{i}(1)];  
    temp_sum = temp_sum + length(proxInd{i});
    new_prox_end    = [ new_prox_end temp_sum];
    len_of_vectors  = [ len_of_vectors length(proxInd{i})];
    
    % subtract one, because arrays in C are zero - based
    new_prox_beg(i) = new_prox_beg(i) - 1;
    new_prox_end(i) = new_prox_end(i) - 1;
end

% substract one from each element - in language C arrays are zero based 
data_to_export.new_prox_beg   = new_prox_beg;
data_to_export.new_prox_end   = new_prox_end;
data_to_export.len_of_vectors = len_of_vectors;

% matrices occuring in the problem
data_to_export.A             = full(dat.A);
data_to_export.Q             = full(dat.Q);
data_to_export.L             = full(L);
data_to_export.matrix        = full(matrix);
data_to_export.scaledL       = full(scaled_L);

% exploiting sparsity: export the indices of the nonzero elements
data_to_export.AI       = indexes_of_nnz_elements(data_to_export.A);
data_to_export.QI       = indexes_of_nnz_elements(data_to_export.Q);
data_to_export.LI       = indexes_of_nnz_elements(data_to_export.L);
data_to_export.matrixI  = indexes_of_nnz_elements(data_to_export.matrix);
data_to_export.scaledLI = indexes_of_nnz_elements(data_to_export.scaledL);

% sizes of matrices containing indices of nonzero elements
data_to_export.AI_ROWS         = size(data_to_export.AI, 1);
data_to_export.AI_COLS         = size(data_to_export.AI, 2);
data_to_export.QI_ROWS         = size(data_to_export.QI, 1);
data_to_export.QI_COLS         = size(data_to_export.QI, 2);
data_to_export.LI_ROWS         = size(data_to_export.LI, 1);
data_to_export.LI_COLS         = size(data_to_export.LI, 2);
data_to_export.matrixI_ROWS    = size(data_to_export.matrixI, 1);
data_to_export.matrixI_COLS    = size(data_to_export.matrixI, 2);
data_to_export.scaledLI_ROWS   = size(data_to_export.scaledLI, 1);
data_to_export.scaledLI_COLS   = size(data_to_export.scaledLI, 2);

% prox types, norm types, constants for projections
data_to_export.prox_type = {};

% process prox.type
for i = 1:nProx
    data_to_export.prox_type{i} = char(prox(i).type);
end

% further processing of the structure prox in order to obtain
% the following data:
%     - normType:   prox(i).dat.p - normBall, normProx
%     - constants:  prox(i).dat.c - ellipse, normball

data_to_export.weight    = [];
data_to_export.constants = [];
data_to_export.normType  = {}; % must be a cell array because of Inf norm
for i = 1:length(data_to_export.prox_type)
    
    % process weights
    data_to_export.weight(i) = full(proxWeight(i));
    
    % check norms
    if( strcmp(data_to_export.prox_type{i}, 'normProx') | strcmp(data_to_export.prox_type{i}, 'normBall') ) 
        data_to_export.normType{i} = num2str(prox(i).dat.p);
    else
        data_to_export.normType{i} = num2str(0);
    end
    
    % check data
    if( strcmp(data_to_export.prox_type{i}, 'ellipse') | strcmp(data_to_export.prox_type{i}, 'normBall') ) 
        data_to_export.constants(i) = prox(i).dat.c;
    else
        data_to_export.constants(i) = 0;
    end
    
end

end

% -----------------------------------------------------------------------
%  Validate input arguments
% -----------------------------------------------------------------------
function validate_inputs(prob, settings)
dat = prob.dat; prox = prob.prox;
assert(all(isfield(dat, {'Q', 'c', 'A', 'b'})) && all(isfield(prox, {'L','l'})), ...
  'split:argchk', 'Problem input data not complete')
assert(isvector(dat.c) && (isvector(dat.b) || isempty(dat.b)), 'split:argchk', 'c and b must be vectors')
dat.c = dat.c(:); dat.b = dat.b(:);
n = size(dat.Q,1);
m = size(dat.A,1);
assert(size(dat.Q,2) == n && size(dat.A,2) == n && size(dat.A,1) == length(dat.b) && length(dat.f) == n, ...
  'split:argchk', 'Problem data matrices have inconsistent sizes')
for i = 1:length(prox)
  assert(isvector(prox(i).l), 'split:argchk', 'prox.l must be a vector')
  prox(i).l = prox(i).l(:);
  assert(size(prox(i).L,2) == n && length(prox(i).l) == size(prox(i).L,1), ...
    'split:argchk', 'Problem data matrices have inconsisten sizes in %i prox operator', i)
end
end