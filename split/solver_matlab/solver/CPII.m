function [sol, stats_out, stats_plot, data_to_export] = CPII(prob, settings)
%
% Implementation of CPII to test the export from the split user
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
if isfield(prob, 'param') && length(prob.param) > 0
  prob = setProbParameters(prob);
end

%%% NOTE: setting relaxation to something other than 1 seems to make SDPs fail

if ~isfield(settings, 'terminationCheckFreq'), settings.terminationCheckFreq = 10; end
if ~isfield(settings, 'relaxation'),           settings.relaxation = 1; end
if ~isfield(settings, 'dualTol'),              settings.dualTol = 1e-3; end
if ~isfield(settings, 'primalTol'),            settings.primalTol = 1e-3; end
if ~isfield(settings, 'maxItr'),               settings.maxItr = 1e3; end
if ~isfield(settings, 'fun_eval'),             settings.fun_eval = 0; end

% Convergence statistics requested
RECORD_STATS = false;
if nargout > 1 
  RECORD_STATS = true;
end

dat   = prob.dat;
prox  = prob.prox;
nProx = length(prox);

% Validate inputs
validate_inputs(prob, settings);

% -----------------------
% Pre-processing
% -----------------------
n = size(dat.Q, 1);
m = size(dat.A, 1);

L = []; l = []; 
nProxVars = zeros(nProx,1);
kk = 1;
for i = 1:nProx
  l            = [l;prox(i).l];
  L            = [L;prox(i).L];
  proxInd{i}   = kk + [0:size(prox(i).L,1)-1];
  kk           = kk + size(prox(i).L,1);
  nProxVars(i) = size(prox(i).L,1);
end
proxWeight = [prob.prox.weight];

% define matrices specific for chambolle-pock
K = [L; dat.A];
matrix_inverse = []; % CANNOT PRECOMPUTE INVERSE

% Bound on K operator
Mu = max(svds(K)); 

% stepsizes
theta = settings.relaxation;
if ~isfield(settings, 'tau'),           settings.tau = 1/Mu; end
tau = settings.tau;
if ~isfield(settings, 'sigma'),         settings.sigma = 1/(tau*Mu^2); end
sigma = settings.sigma; % tau and sigma set according to Theorem 2 of paper
if ~isfield(settings, 'gamma'),         settings.gamma = 0.1*min(eig(dat.Q)); end
gamma = settings.gamma; % Empirically, the algorithm behaves better for gamma = 0.1*lambda_min(Q)

% % for i = 1:nProx
% %     % special cases: normBall, ellipse - I need to scale the proxWeight vector
% %     % Perhaps it is not the best solution, but then I need to change the 
% %     % interface, i.e. add one more field to prob structure e.g. sigma,
% %     % however, it is a specific case occuring only in chambolle
% %     if strcmp(char(prox(i).typeConj), 'normProx')
% %          proxWeight(i) = sigma * prox(i).dat.c;
% %     end
% %     if strcmp(char(prox(i).typeConj),  'ellipseConj')
% %          proxWeight(i) = sigma * sqrt(prox(i).dat.c); 
% %     end
% %         
% % end

%  Variables initialization and warm-starting + export
if ~isfield(settings, 'primal_vars'),     settings.primal_vars = zeros(n, 1); end
x = settings.primal_vars; 
x_bar = settings.primal_vars;
data_to_export.x = x; data_to_export.xbar = x_bar;

if ~isfield(settings, 'dual_vars_nu'),    settings.dual_vars_nu = zeros(size(dat.A, 1), 1); end
nu = settings.dual_vars_nu;
data_to_export.nu = nu;

if ~isfield(settings, 'dual_vars_p'),     settings.dual_vars_p = zeros(size(L, 1), 1); end
p = settings.dual_vars_p;    
data_to_export.p = p;
p_nu = [p; nu];

% Scaling of the residuals
if ~isfield(settings, 'scale_rPrimal'),   settings.scale_rPrimal = eye(n); end
SP = settings.scale_rPrimal;
if ~isfield(settings, 'scale_rDual'),    settings.scale_rDual = eye(size(K,1)); end
SD = settings.scale_rDual;
invSP = inv(SP);
invSD = inv(SD);

% export them here: because they will be changed, gamma not though - I
% think
data_to_export.theta             = theta;
data_to_export.sigma             = sigma;
data_to_export.tau               = tau;
data_to_export.gamma             = gamma;

% Main
iii = 0;
obj = 0;
for jjj = 1:settings.maxItr
  iii = iii + 1;
    
  % step 0
  p_prev    = p;
  nu_prev   = nu;
  x_prev    = x;
  x_bar_prev = x_bar;
  p_nu_prev = [p_prev; nu_prev];
  
  % step 1: compute p and nu
  p = p + sigma*(L*x_bar + l);
  nu = nu + sigma*(dat.A*x_bar - dat.b);
     
  % step 2: compute projection
  for i = 1:nProx
    % special cases: normBall, ellipse - I need to scale the proxWeight vector
    % Perhaps it is not the best solution, but then I need to change the 
    % interface, i.e. add one more field to prob structure e.g. sigma,
    % however, it is a specific case occuring only in chambolle
    if strcmp(char(prox(i).typeConj), 'normProx')
         proxWeight(i) = sigma * prox(i).dat.c;
    end
    if strcmp(char(prox(i).typeConj),  'ellipseConj')
         proxWeight(i) = sigma * sqrt(prox(i).dat.c); 
    end
      ind = proxInd{i};
      p(ind) = prox(i).conjFunc(p(ind), proxWeight(i)); 
  end
  p_nu = [p; nu];
  
  % step 3: compute x
  % this formula represents the conventional procedure
  % x = (dat.Q + eye(size(dat.Q, 1))*(1/tau)) \ ((1/tau)*x - dat.f' - (K'*p_nu));
  
  % exploit structure of matrix Q - simpler C implementation
  % pull out elements from the diagonal and update x
  diagQ    = diag(dat.Q);
  warning('Implementation in the current form can only handle the diagonal Q')
  diagItau = diag(eye(size(dat.Q, 1))*(1/tau));
  inverse  = 1./(diagQ + diagItau);
  x = inverse .* ((1/tau)*x - dat.f' - (K'*p_nu));
  
  % step 4: update stepsize
  theta = 1 / (sqrt(1+2*gamma*tau));
  tau = theta*tau;
  sigma = sigma / theta;
      
  % step 5: update x_bar
  x_bar = x + theta*(x - x_prev);
  
  % function evaluation
  if settings.fun_eval == 1
      obj = .5*(invSP*x)'*dat.Q*(invSP*x) + dat.f*(invSP*x);
  end
  
  % step 6: convergence check
%   scaledK = theta*-K;  
  r  = (1/tau)*(invSP*(x_prev-x)) - (K*SP)'*(p_nu_prev - p_nu); 
  s  = (1/sigma)*SD'*(p_nu_prev-p_nu) - invSD*K*(x - x_bar);

  rPrimal = norm(r);
  rDual   = norm(s);
  
  if RECORD_STATS
     if settings.fun_eval == 1
        stats_plot(iii).f_eval  = obj;
     end
    stats_plot(iii).rDual       = rDual;
    stats_plot(iii).rPrimal     = rPrimal;
    stats_plot(iii).p           = p;
    stats_plot(iii).x           = x;
    stats_plot(iii).nu          = nu;
    stats_plot(iii).s           = s;
    stats_plot(iii).r           = r;
    stats_plot(iii).tau         = tau;
    stats_plot(iii).sigma       = sigma;
    stats_plot(iii).theta       = theta;
  end 
  
  if rDual < settings.dualTol && rPrimal < settings.primalTol
    cprintf([0.1 0.5 0.1], '\n\n>>>>> Stopping on optimality after %i steps <<<<<\n\n', iii);
    break
  end

end

% record outside the loop hence we have only the most recent values
if RECORD_STATS
   stats_out.rDual   = rDual;
   stats_out.rPrimal = rPrimal;
   stats_out.p       = p;
   stats_out.x       = x;
   stats_out.nu      = nu;
   stats_out.s       = s;
   stats_out.r       = r;
   stats_out.numiter = jjj;
end

sol.problem = 0;
if rDual > settings.dualTol || rPrimal > settings.primalTol
  cprintf('err', 'Failed to reach optimality after %i iterations\n', settings.maxItr);
  sol.problem = 1;
end

sol.x  = x;
sol.p  = p;
sol.nu = nu;

% prepare data for export to C

% common data required to solve the problem, e.g. stepsize, number of
% states, inputs ...
% data_to_export.num_of_opt_var    = n;
% data_to_export.num_of_states     = m;
% data_to_export.nProx             = nProx;
% data_to_export.settings          = settings;
% 
% % sizes of vectors
% data_to_export.b_LEN          = length(dat.b);
% data_to_export.f_LEN          = length(dat.f);
% data_to_export.c_LEN          = length(dat.c);
% data_to_export.l_LEN          = length(l);
% data_to_export.x_LEN          = length(x);
% data_to_export.nu_LEN         = length(nu);
% data_to_export.p_LEN          = length(p);
% data_to_export.p_nu_LEN       = length(p_nu);
% data_to_export.diagQ_LEN      = length(diag(dat.Q));
% 
% % sizes of matrices
% data_to_export.A_ROWS         = size(dat.A, 1);
% data_to_export.A_COLS         = size(dat.A, 2);
% data_to_export.L_ROWS         = size(L, 1);
% data_to_export.L_COLS         = size(L, 2);
% data_to_export.K_ROWS         = size(K, 1);
% data_to_export.K_COLS         = size(K, 2);
% data_to_export.KT_ROWS        = size(K', 1);
% data_to_export.KT_COLS        = size(K', 2);
% data_to_export.SP_ROWS        = size(SP, 1);
% data_to_export.SP_COLS        = size(SP, 2);
% data_to_export.SD1_ROWS       = size(SD, 1);
% data_to_export.SD1_COLS       = size(SD, 2);
% 
% % vectors
% % data_to_export.x              = zeros(n,1);
% data_to_export.x_prev         = zeros(n,1);
% % data_to_export.xbar           = zeros(n,1);
% data_to_export.xbar_prev      = zeros(n,1);
% % data_to_export.nu             = zeros(size(dat.A, 1),1);
% data_to_export.nu_prev        = zeros(size(dat.A, 1),1);
% % data_to_export.p              = zeros(size(L, 1),1);
% data_to_export.p_prev         = zeros(size(L, 1),1);
% data_to_export.p_nu           = zeros(length(p_nu), 1);
% data_to_export.p_nu_prev      = zeros(length(p_nu), 1);
% data_to_export.b              = full(dat.b);
% data_to_export.f              = full(dat.f);
% data_to_export.c              = full(dat.c);
% data_to_export.l              = full(l);
% data_to_export.diagQ          = full(diag(dat.Q));
% 
% new_prox_beg   = [];
% new_prox_end   = [];
% len_of_vectors = [];
% 
% temp_sum = 0;
% for i = 1:nProx
%     new_prox_beg    = [ new_prox_beg proxInd{i}(1)];  
%     temp_sum = temp_sum + length(proxInd{i});
%     new_prox_end    = [ new_prox_end temp_sum];
%     len_of_vectors  = [ len_of_vectors length(proxInd{i})];
%     
%     % subtract one from each element, because arrays in C are zero - based
%     new_prox_beg(i) = new_prox_beg(i) - 1;
%     new_prox_end(i) = new_prox_end(i) - 1;
% end
% 
% % substract one from each element - in language C arrays are zero based 
% data_to_export.new_prox_beg   = new_prox_beg;
% data_to_export.new_prox_end   = new_prox_end;
% data_to_export.len_of_vectors = len_of_vectors;
% 
% % matrices occuring in the problem
% data_to_export.A             = full(dat.A);
% data_to_export.K             = full(K);
% data_to_export.KT            = full(K');
% data_to_export.L             = full(L);
% data_to_export.SP            = SP;
% data_to_export.SD1           = SD;
% 
% % exploiting sparsity: export the indices of the nonzero elements
% data_to_export.AI        = indexes_of_nnz_elements(data_to_export.A);
% data_to_export.KI        = indexes_of_nnz_elements(data_to_export.K);
% data_to_export.KTI       = indexes_of_nnz_elements(data_to_export.K');
% data_to_export.LI        = indexes_of_nnz_elements(data_to_export.L);
% data_to_export.SPI       = indexes_of_nnz_elements(SP);
% data_to_export.SD1I      = indexes_of_nnz_elements(SD);
% 
% % sizes of matrices containing indices of nonzero elements
% data_to_export.AI_ROWS         = size(data_to_export.AI, 1);
% data_to_export.AI_COLS         = size(data_to_export.AI, 2);
% data_to_export.KI_ROWS         = size(data_to_export.KI, 1);
% data_to_export.KI_COLS         = size(data_to_export.KI, 2);
% data_to_export.KTI_ROWS        = size(data_to_export.KTI, 1);
% data_to_export.KTI_COLS        = size(data_to_export.KTI, 2);
% data_to_export.LI_ROWS         = size(data_to_export.LI, 1);
% data_to_export.LI_COLS         = size(data_to_export.LI, 2);
% data_to_export.SPI_ROWS        = size(data_to_export.SPI , 1);
% data_to_export.SPI_COLS        = size(data_to_export.SPI , 2);
% data_to_export.SD1I_ROWS       = size(data_to_export.SD1I, 1);
% data_to_export.SD1I_COLS       = size(data_to_export.SD1I, 2);
% 
% % prox types, norm types, constants for projections
% % Since it is chambolle we need the conjugates instead of the original prox
% % operators
% data_to_export.proxConj_type = {};
% 
% % process prox.type
% for i = 1:nProx
%     data_to_export.proxConj_type{i} = char(prob.prox(i).typeConj);
% end
% 
% % further processing of the structure prox in order to obtain
% % the following data:
% %     - normType:   prox(i).dat.pDual - normBall, normProx, ellipseConj
% %     - constants:  prox(i).dat.c     - ellipseConj, normball
% % 
% data_to_export.weight    = [];
% data_to_export.constants = [];
% data_to_export.dualNormType  = {}; % must be a cell array because of Inf norm
% for i = 1:length(data_to_export.proxConj_type)
%     
%     % process weights
%     data_to_export.weight(i) = full(proxWeight(i));
%     
%     % check norms
%     if( strcmp(data_to_export.proxConj_type{i}, 'normProx') | strcmp(data_to_export.proxConj_type{i}, 'normBall') | strcmp(data_to_export.proxConj_type{i}, 'ellipseConj') ) 
%         data_to_export.dualNormType{i} = num2str(prox(i).dat.pDual);
%     else
%         data_to_export.dualNormType{i} = num2str(0);
%     end
%     
%     % check data
%     if( strcmp(data_to_export.proxConj_type{i}, 'ellipseConj') | strcmp(data_to_export.proxConj_type{i}, 'normProx') | strcmp(data_to_export.proxConj_type{i}, 'normBall') ) 
%         data_to_export.constants(i) = prox(i).dat.c;
%     else
%         data_to_export.constants(i) = 0;
%     end
%     
%     % special case: norm in the objective:
%     % if 
%     
% end

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