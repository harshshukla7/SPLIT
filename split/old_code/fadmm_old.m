function [sol, stats_out, stats_plot] = fadmm(prob, settings)
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
if ~isfield(settings, 'dualTol'),              settings.dualTol = 1e-3; end
if ~isfield(settings, 'primalTol'),            settings.primalTol = 1e-3; end
if ~isfield(settings, 'maxItr'),               settings.maxItr = 1e3; end
if ~isfield(settings, 'fun_eval'),             settings.fun_eval = 0; end

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
rho   = settings.rho;

% -----------------------
% Pre-processing
% -----------------------

n = size(dat.Q,1);
m = size(dat.A,1);

% Setup the KKT matrix and pre-solve
Q = dat.Q;
sum_l = zeros(n,1);
sum_L = [];
L = []; l = []; 
nProxVars = zeros(nProx,1);
kk = 1;
for i = 1:nProx
  Q            = Q + rho*prox(i).L'*prox(i).L;
  sum_l        = sum_l + rho*prox(i).L'*prox(i).l;
  sum_L        = [sum_L prox(i).L'];
  l            = [l;prox(i).l];
  L            = [L;prox(i).L];
  proxInd{i}   = kk + [0:size(prox(i).L,1)-1];
  kk           = kk + size(prox(i).L,1);
  nProxVars(i) = size(prox(i).L,1);
end
proxWeight = [prob.prox.weight]; % Weights for the prox functions

% Solution to KKT system
if (sum(any(Q-diag(diag(Q))))==0) && (isempty(dat.A))
    Lkkt = eye(n+m); Pkkt = eye(n+m); Dkkt = Q;
else
   KKT = [Q dat.A'; dat.A zeros(m)];
   [Lkkt,Dkkt,Pkkt] = ldl(KKT);
end
rhsKKTconst = [-dat.f' - sum_l; dat.b];
rhsKKTvar   = [-rho*sum_L;sparse(length(dat.b),size(sum_L,2))];

% Form the augmented lagrangian
% min 0.5*x'*Q*x + f'*x + sum_i Fi(y_i) + rho/2*sum_i ||Li*x + li - y_i + lam_i||
%      A*x = b

% % % Matrices to store the solution
% % x       = zeros(n,1);
% % y       = zeros(sum(nProxVars),1);
hat.y   = zeros(sum(nProxVars),1); % hat quantities refer to accelerated sequencies
% % lam     = zeros(sum(nProxVars),1);
hat.lam = zeros(sum(nProxVars),1);
% % q       = zeros(size(L,1),1);
rDual   = 0;
rPrimal = 0;
alpha   = 1;
eta     = 0.999; % eta~1 fewer restarts

%  Variables initialization and warm-starting
if ~isfield(settings, 'primal_vars_x'),     settings.primal_vars_x = zeros(n, 1); end
x = settings.primal_vars_x; 
if ~isfield(settings, 'primal_vars_y'),    settings.primal_vars_y = zeros(size(L,1),1); end
y = settings.primal_vars_y;
if ~isfield(settings, 'dual_vars_lambda'),     settings.dual_vars_lambda = zeros(size(L,1),1); end
lam = settings.dual_vars_lambda;         

% Scaling of the residuals
if ~isfield(settings, 'scale_rPrimal'),   settings.scale_rPrimal = eye(n); end
SP = settings.scale_rPrimal;
invSP = inv(SP);
if ~isfield(settings, 'scale_rDual'),   settings.scale_rDual = eye(size(L,1)); end
SD = settings.scale_rDual;
invSD = inv(SD);

% precompute some common expression
% scaled_L = -rho*sum_L;
const_vector = Pkkt*(Lkkt'\(Dkkt\(Lkkt\(Pkkt'* rhsKKTconst))));
matrix = Pkkt*(Lkkt'\(Dkkt\(Lkkt\(Pkkt' * rhsKKTvar))));

c = norm(L*x+l-y)^2;
iii = 0;
for jjj = 1:settings.maxItr
  iii = iii + 1;

  % vector needed to check convergence
  prev.y = y;

  % Step 1 : solve linear system
  t = const_vector + matrix * (hat.lam/rho - hat.y);
  x = t(1:n);
  
  % Step 2 : compute prox operators
  q = L*x + l + 1/rho*hat.lam;
 
  for i = 1:nProx
    ind = proxInd{i};
    y(ind) = prox(i).func(q(ind), 1/rho*proxWeight(i));
  end

  % Update the lagrange multiplier
  prev.lam = lam;
  lam = hat.lam + rho*(L*x + l - y);
  
  % function evaluation
  if settings.fun_eval == 1
      obj = .5*(invSP*x)'*dat.Q*(invSP*x) + dat.f*(invSP*x);
  end
   
  % Convergence checks
  s = rho*L'*(prev.y-y);
  r = (L*x + l - y);
  
  rDual   = norm(s);
  rPrimal = norm(r);
  
  % Adaptive restart of Nesterov's relaxation
  prev.c = c;
  c = 1/rho * norm(lam-hat.lam)^2 + rho*norm(y-hat.y)^2;
  
  if c < eta*prev.c
      prev.alpha = alpha;
      alpha = (1+sqrt(4*prev.alpha^2+1))/2;    
      hat.y = y + (prev.alpha-1)*(y - prev.y)/alpha;
      hat.lam = lam + (prev.alpha-1)*(lam - prev.lam)/alpha;
  else
      alpha = 1; 
      hat.lam = lam;
      hat.y = y;
      c = 1/eta * prev.c;
  end
  
    % we need record stats inside the loop as well for plots
  if RECORD_STATS
     if settings.fun_eval == 1
        stats_plot(iii).f_eval  = obj;
     end
     stats_plot(iii).rDual   = rDual;
     stats_plot(iii).rPrimal = rPrimal;
     stats_plot(iii).x       = x;
     stats_plot(iii).y       = y;
     stats_plot(iii).lam     = lam;
     stats_plot(iii).s       = s;
     stats_plot(iii).r       = r;
     stats_plot(iii).rho     = rho;
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
   stats_out.x       = x;
   stats_out.y       = y;
   stats_out.lam     = lam;
   stats_out.s       = s;
   stats_out.r       = r;
   stats_out.rho     = rho;
   stats_out.numiter = jjj;
end

sol.problem = 0;
if rDual > settings.dualTol || rPrimal > settings.primalTol
  cprintf('err', 'Failed to reach optimality after %i iterations\n', settings.maxItr);
  sol.problem = 1;
end

sol.lam = lam;
sol.x = x;
sol.y = y;

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
