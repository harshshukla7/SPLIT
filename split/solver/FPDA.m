function [sol, stats_out, stats_plot] = FPDA(prob, settings)
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

%% Scalings, relaxations and stepsizes
if ~isfield(settings, 'precondition') % No prec, No adapt
    settings.scale_rPrimal = eye(n);
    settings.scale_rDual = eye(size(L,1));
    Mu = max(svds(L))^2; % less conservative bound
% % %     if ~isempty(dat.Q) % existence of smooth term changes stepsize
% % %         settings.tau = 0.05/max(eig(dat.Q));
% % %         tau = settings.tau;
% % %         settings.rho = (1-tau*max(eig(dat.Q))/2)^2 / (tau*Mu);
% % %         rho = settings.rho;
% % %     else
% % %         settings.tau = 1/sqrt(Mu);
% % %         tau = settings.tau;
% % %         settings.rho = tau;
% % %         rho = settings.rho;
% % %     end
    
end
prec.Dinv = settings.scale_rPrimal; % D^(-1) (preconditioner)
prec.D = inv(settings.scale_rPrimal); % D
prec.T = settings.scale_rPrimal*settings.scale_rPrimal'; % T = D^(-1)*D^(-T)
prec.P  = settings.scale_rDual'*settings.scale_rDual; % P=E'E
prec.E  = settings.scale_rDual; % E (preconditioner)
prec.L  = prec.E*L*prec.Dinv;
% % if isfield(settings, 'precondition') 
% %     Mu      = settings.Mu;
% %    
% % end

%% stepsizes
if ~isfield(settings, 'gamma'),         settings.gamma = 0.05*min(eig(prec.Dinv'*dat.Q*prec.Dinv)); end
gamma = settings.gamma; % strong convexity constant
if ~isfield(settings, 'Lip'),           settings.Lip = max(eig(prec.Dinv'*dat.Q*prec.Dinv)); end
Lip = settings.gamma; % Lipschitz constant for the gradient
if ~isfield(settings, 'kappa'),         settings.kappa = Lip+1; end
kappa = settings.kappa; 
if ~isfield(settings, 'tau'),           settings.tau = 0.04*(2*gamma/Lip); end
tau = settings.tau;
if ~isfield(settings, 'theta'),         settings.relaxation = 1/sqrt(1+tau*(2*gamma-Lip*tau)/kappa); end
theta = settings.relaxation;
if ~isfield(settings, 'rho'),           settings.rho = 1 / (tau*theta*Mu); end
rho = settings.rho;
    

%%  Variables initialization and warm-starting + export
if ~isfield(settings, 'primal_vars'),     settings.primal_vars = zeros(n, 1); end
x = settings.primal_vars;
p = x;

if ~isfield(settings, 'dual_vars_lam'),     settings.dual_vars_lam = zeros(size(L, 1), 1); end
lam = settings.dual_vars_lam; 
lamd = prec.E*lam;

%% Main
iii = 0;
obj = 0;
for jjj = 1 : settings.maxItr
  iii = iii + 1;
 
  %% step 1: linear system solve
  y = (eye(size(dat.Q))-(tau/kappa)*prec.T*dat.Q)*p-(tau/kappa)*prec.T*(L'*lamd+dat.f');
  prev.p    = p; 
  if ~isempty(dat.A) % dynamics update
      p = y + dat.A'*((dat.A*dat.A')\(dat.b-dat.A*y));
  else
      p = y;
  end
  prev.x = x;
  x = p + theta*(p-prev.p);

  %% step 2: compute projection
  prev.lamd = lamd;
  lamd_tilde = lamd + rho*prec.P*(L*x + l);
  for i = 1:nProx
    if strcmp(char(prox(i).typeConj), 'normProx')
%          proxWeight(i) = rho*prec.P( proxInd{i}, proxInd{i}) * prox(i).dat.c;
         proxWeight(i) = rho * prox(i).dat.c;
    end
    if strcmp(char(prox(i).typeConj),  'ellipseConj')
%          proxWeight(i) = rho*prec.P( proxInd{i}, proxInd{i}) * sqrt(prox(i).dat.c); 
         proxWeight(i) = rho * sqrt(prox(i).dat.c); 
    end
      ind = proxInd{i};
      lamd(ind) = prox(i).conjFunc(lamd_tilde(ind), proxWeight(i));
  end

% %   % step 6: convergence check
% %   % ATTENTION! data unscaled, x unscaled, p,mu scaled 
  dx = x-prev.x;
  dlamd = lamd-prev.lamd;
  r = (dat.Q-(1/tau)*prec.D)*dx + prec.Dinv'*(L'*(dlamd));
  s = (1/rho)*(prec.E'\(-dlamd)) + prec.E*L*dx;
  
%   L*dx
%   rho
%   1/rho
  rPrimal = norm(r,1);
  rDual   = norm(s,1);

  % step 4: update stepsize
  tau = theta*tau;
  theta = 1/sqrt(1+tau*(2*gamma-Lip*tau)/kappa);
  rho = rho / theta;
  
  if RECORD_STATS
     if settings.fun_eval == 1
        stats_plot(iii).f_eval  = obj;
     end
    stats_plot(iii).rDual       = rDual;
    stats_plot(iii).rPrimal     = rPrimal;
    stats_plot(iii).lam         = prec.E\lamd;
    stats_plot(iii).x           = x;
    stats_plot(iii).s           = s;
    stats_plot(iii).r           = r;
    stats_plot(iii).tau         = tau;
    stats_plot(iii).rho         = rho;
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
   stats_out.lam     = prec.E\lamd;
   stats_out.x       = x;
   stats_out.s       = s;
   stats_out.r       = r;
   stats_out.tau     = tau;
   stats_out.rho     = rho;
   stats_out.numiter = jjj;
end

sol.problem = 0;
if rDual > settings.dualTol || rPrimal > settings.primalTol
  cprintf('err', 'Failed to reach optimality after %i iterations\n', settings.maxItr);
  sol.problem = 1;
end

sol.x    = x;
sol.lam  = prec.E\lamd;
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