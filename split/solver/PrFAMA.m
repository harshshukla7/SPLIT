function [sol, stats_out, stats_plot] = PrFAMA(prob, settings)
%
% Quick implementation of AMA with optional scaling, Nesterov acceleration
% and restart (not supported yet)
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
% If the problem is parametric, then the current values of the parameters
% are used when solving the problem. Use set to set the parameter values
%
% Structure 'prec' containts the scaling matrices E, D, P, T that should
% always be positive definite diagonals. If not specified set to identity
%
% Variables followed by a 'd' indicate scaling by constraint preconditioner E
%

if nargin < 2
  settings = [];    
end

% The problem is parametric
if isfield(prob, 'param') && length(prob.param) > 0
  prob = setProbParameters(prob);
end

if ~isfield(settings, 'terminationCheckFreq'), settings.terminationCheckFreq = 10; end
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

%% Data
n = size(dat.Q,1);
m = size(dat.A,1);

% Extract Lx + l in matrix form
L = []; l = [];
for i = 1:nProx
  l            = [l;prox(i).l];
  L            = [L;prox(i).L];
end

%% Scalings
if ~isfield(settings, 'scale_rPrimal')
    settings.scale_rPrimal = eye(n); 
end
prec.Dinv = settings.scale_rPrimal;
prec.T = settings.scale_rPrimal'*settings.scale_rPrimal;
if ~isfield(settings, 'scale_rDual') 
    settings.scale_rDual = eye(size(L,1)); 
end
prec.P  = settings.scale_rDual'*settings.scale_rDual; % P=E'E
prec.E  = settings.scale_rDual; % E (preconditioner)

% % % Primal change of basis
% % if ~isfield(settings, 'changeBasisPrimal')
% %     settings.changeBasisPrimal = eye(n); 
% % end
% % SP = settings.changeBasisPrimal;

%% Setup the KKT matrix and pre-solve
Q = dat.Q;
sum_L = [];
L = []; l = []; 
nProxVars = zeros(nProx,1);
kk = 1;

for i = 1:nProx
  proxInd{i}   = kk + [0:size(prox(i).L,1)-1];
  ind          = proxInd{i};
  sum_L        = [sum_L prox(i).L'*prec.E(ind,ind)'];
  l            = [l;prec.E(ind,ind)*prox(i).l]; % ATTENTION! L and l are left scaled
  L            = [L;prec.E(ind,ind)*prox(i).L];
  proxInd{i}   = kk + [0:size(prox(i).L,1)-1];
  kk           = kk + size(prox(i).L,1);
  nProxVars(i) = size(prox(i).L,1);
end
proxWeight = [prob.prox.weight]; % Weights for the prox functions

% Solution to KKT system
KKT = [Q dat.A'; dat.A zeros(m)];
[Lkkt,Dkkt,Pkkt] = ldl(KKT,1e-6);
rhsKKTconst = [-dat.f'; dat.b];
rhsKKTvar   = [sum_L;sparse(length(dat.b),size(sum_L,2))];

%% Compute the largest stepsize
if isempty(dat.A)
    sigma_f = min(eig(prec.Dinv'*dat.Q*prec.Dinv));
    if sigma_f <= 1e-8; % tolerance
        error('Cannot use FAMA - objective not strongly convex')
    end
else
    Nul = null(full(dat.A));
    if min(eig(Nul'*prob.dat.Q*Nul)) <= 1e-8;
        error('Cannot use FAMA - objective not strongly convex in the nullspace of the dynamics')
    end
    sigma_f = min(eig(Nul'*prec.Dinv'*dat.Q*prec.Dinv*Nul));
end
eigenvalue_LL = max(eig(L'*L)); % |L|_2^2

if ~isfield(settings, 'rho'),     settings.rho = sigma_f/eigenvalue_LL; end
rho = settings.rho;

%%  Variables initialization and warm-starting
if ~isfield(settings, 'primal_vars_x'),     settings.primal_vars_x = zeros(n, 1); end
x = settings.primal_vars_x; 

if ~isfield(settings, 'primal_vars_y'),    settings.primal_vars_y = zeros(size(L,1),1); end
yd = prec.E*settings.primal_vars_y;

if ~isfield(settings, 'dual_vars_lam'),     settings.dual_vars_lam = zeros(size(L,1),1); end
lam = settings.dual_vars_lam;         
prev.lam = lam;

% -----------------------
% Main
% -----------------------
beta = 1;
prev.beta = beta;
iii = 0;
for jjj = 1  :  settings.maxItr
  iii = iii + 1;
  
  %% Step 1: Nesterov acceleration
  beta = (1+sqrt(4*prev.beta^2+1))/2;
  prev.beta = beta;
  
  
  hat.lam = lam + (prev.beta-1)*(lam - prev.lam)/beta;
  
%   disp('lambda hat after nesterov')
%   hat.lam
  
  
  %% Step 2 : solve linear system
  prev.x = x;
  
%   disp('before kkt solve')
%   (rhsKKTconst - rhsKKTvar*hat.lam)
  
  t = KKT \ (rhsKKTconst - rhsKKTvar*hat.lam);
  
%   disp('kkt solve')
%   t
  
  x = t(1:n);
  
  %% Step 3 : compute prox operators
  q =  prec.E\(L*x + l + (1/rho*hat.lam));

  
  for i = 1:nProx
    ind = proxInd{i};
    yd(ind) = diag(prec.E(ind,ind)).*prox(i).func(q(ind), 1/rho*proxWeight(i)); 
  end
  
%   disp('prox solve')
%   yd     
  
  %rho
  %L*x+l+hat.lam/rho
  %% Step 4: update the Lagrange multiplier
  prev.lam = lam;
  
%   disp('rho is')
%           rho
%           disp('workdual is')
%           (L*x+l) + hat.lam/rho
          
%   disp('lam solve')    
  lam = hat.lam + rho*(L*x + l - yd);
  
  %% Step 5: adaptive restart
  if isfield(settings,'restart')
      % generalized gradient evaluation
      test = (hat.lam-lam)' * (lam-prev.lam);
      if test > 0 % monotonicity test
          hat.lam = prev.lam;
          
          %pause
          t = KKT \ (rhsKKTconst - rhsKKTvar*hat.lam);
          x = t(1:n);
          %pause
          q =  prec.E\(L*x + l + (1/rho*hat.lam));
          for i = 1:nProx
            ind = proxInd{i};
            yd(ind) = diag(prec.E(ind,ind)).*prox(i).func(q(ind), 1/rho*proxWeight(i)); 
          end
          
            
          %pause
          prev.lam = hat.lam;
          lam = hat.lam + rho*(L*x + l - yd);
          %pause
      end
  end
 
  
% %   %% Convergence checks on original
% %   % L, l scaled, z, lam unscaled, yd scaled
% %   s = (prec.E\L)'*(prev.lam - lam);
% %   r = prec.E\(L*x + l - yd);

  %% step 4: convergence check
  % L, l scaled, z, lam unscaled, yd scaled
  s = (L*prec.Dinv)'*(prev.lam - lam);
  r = L*x + l - yd;
  
  rPrimal = norm(r,1);
  rDual   = norm(s,1);
  
   
  %% record statistics (optional)
  if RECORD_STATS
     if settings.fun_eval == 1
        stats_plot(iii).f_eval  = obj;
     end
     stats_plot(iii).rDual   = rDual;
     stats_plot(iii).rPrimal = rPrimal;
     stats_plot(iii).x       = x;
     stats_plot(iii).y       = prec.E\yd;
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
   stats_out.y       = prec.E\yd;
   stats_out.lam     = lam;
   stats_out.s       = s;
   stats_out.r       = r;
   stats_out.rho     = rho;
   stats_out.numiter = jjj;
   stats_out.KKT     = KKT;
end

sol.problem = 0;
if rDual > settings.dualTol% || rPrimal > settings.primalTol
  cprintf('err', 'Failed to reach optimality after %i iterations\n', settings.maxItr);
  sol.problem = 1;
end

sol.lam = lam;
sol.x = x;
sol.y = prec.E\yd;
sol.rho = rho;

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