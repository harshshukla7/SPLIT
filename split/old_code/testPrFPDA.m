function [sol, stats_out, stats_plot] = testPrFPDA(prob, settings)
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
Mu1 = 0;
for i = 1:nProx
  l            = [l;prox(i).l];
  L            = [L;prox(i).L];
  proxInd{i}   = kk + [0:size(prox(i).L,1)-1];
  kk           = kk + size(prox(i).L,1);
  nProxVars(i) = size(prox(i).L,1);
end
proxWeight = [prob.prox.weight];

%% Scalings
if ~isfield(settings, 'scale_rPrimal')
    settings.scale_rPrimal = eye(n); 
% else
%     tau = 1;
end
prec.Dinv = settings.scale_rPrimal; % D^(-1) (preconditioner)
prec.T = settings.scale_rPrimal*settings.scale_rPrimal'; % T = D^(-1)*D^(-T)
if ~isfield(settings, 'scale_rDual') 
    if ~isempty(dat.A)
        settings.scale_rDual = eye(size(L,1)+n); 
    else
        settings.scale_rDual = eye(size(L,1)); 
    end
% else
%     sigma = 1;
end
prec.P = settings.scale_rDual'*settings.scale_rDual; % P=E'E
prec.E = settings.scale_rDual; % E (preconditioner)

% % % if ~isfield(settings, 'changeBasisPrimal')
% % %     settings.changeBasisPrimal = eye(n); 
% % % end
% % % SP = settings.changeBasisPrimal;

% Bound on K operator 
if ~isempty(dat.A)
    prec.E1 = prec.E(1:size(L,1),1:size(L,1));
    prec.E2 = prec.E(size(L,1)+1:size(L,1)+n,size(L,1)+1:size(L,1)+n);
    prec.P1 = prec.P(1:size(L,1),1:size(L,1));
    prec.P2 = prec.P(end-n+1:end, end-n+1:end);
    prec.L  = prec.E1*L*prec.Dinv;
    prec.I  = prec.E2*eye(n)*prec.Dinv;
    prec.K   = [prec.L; prec.I];
    Mu       = max(svds(prec.K))^2; % sigma(prec.K'*prec.K)^2 = |prec.K'*prec.K|
else 
    prec.E1 = prec.E(1:size(L,1),1:size(L,1));
    prec.E2 = zeros(n,n);
    prec.P1 = prec.P(1:size(L,1),1:size(L,1));
    prec.P2 = zeros(n, n);
    prec.L = prec.E1*L*prec.Dinv;
    prec.K = prec.L;
    Mu     = max(svds(prec.K))^2;
end

%% stepsizes
if ~isfield(settings, 'gamma'),         settings.gamma = 1e-1*min(eig(prec.Dinv'*dat.Q*prec.Dinv)); end
gamma = settings.gamma; % Empirically, the algorithm behaves better for gamma = 0.1*lambda_min(Q)
if ~isfield(settings, 'lam'),           settings.lam = 1; end
lam = settings.lam; % Empirically, the algorithm behaves better for gamma = 0.1*lambda_min(Q)
if ~isfield(settings, 'tau'),           settings.tau = gamma; end
tau = settings.tau;
if ~isfield(settings, 'sigma'),         settings.sigma = 1 / (tau*Mu); end
sigma = settings.sigma; % tau and sigma set according to Theorem 2 of paper
% %         settings.tau = 0.95/sqrt(Mu);
% %         tau = settings.tau;
% %         settings.sigma = 0.95/(tau*(Mu));
% %         sigma = settings.sigma;
if ~isfield(settings, 'theta'),         settings.relaxation = 1/sqrt(1+tau*2*gamma/lam); end
theta = settings.relaxation;
    

%%  Variables initialization and warm-starting + export
if ~isfield(settings, 'primal_vars'),     settings.primal_vars = zeros(n, 1); end
xp = settings.primal_vars; 
xp_bar = settings.primal_vars;

if ~isfield(settings, 'dual_vars_lambda'),     settings.dual_vars_lambda = zeros(size(L, 1), 1); end
lambda = settings.dual_vars_lambda; 

if ~isfield(settings, 'dual_vars_mu'),    settings.dual_vars_mu = zeros(size(dat.A, 1), 1); end
mu = settings.dual_vars_mu;

%% Main
iii = 0;
obj = 0;
mu = xp; % left scaled
for jjj = 1:settings.maxItr
  iii = iii + 1;
    
  % step 0
  prev.xp    = xp;
  
  % step 1: linear system solve
  xp = prec.Dinv\...
      ( ((eye(size(dat.Q))+tau*prec.T*dat.Q)) \ (prec.Dinv*prev.xp-tau*prec.T*((prec.E1*L)'*lambda+prec.E2'*mu+dat.f')) );
  xp_bar = 2*xp - prev.xp;
     
  tilde.lambda = prec.E1'*lambda + sigma*prec.P1*(L*prec.Dinv*xp_bar + l);

  % step 2: compute projection
  prev.lambda = lambda;
  for i = 1:nProx
    if strcmp(char(prox(i).typeConj), 'normProx')
%          proxWeight(i) = sigma*prec.P( proxInd{i}, proxInd{i}) * prox(i).dat.c;
         proxWeight(i) = sigma * prox(i).dat.c;
    end
    if strcmp(char(prox(i).typeConj),  'ellipseConj')
%          proxWeight(i) = sigma*prec.P( proxInd{i}, proxInd{i}) * sqrt(prox(i).dat.c); 
         proxWeight(i) = sigma * sqrt(prox(i).dat.c); 
    end
      ind = proxInd{i};
      lambda(ind) = prec.E1\(prox(i).conjFunc(tilde.lambda(ind), proxWeight(i))); 
  end

   if ~isempty(dat.A)
      prev.mu    = mu;
      mu = sigma*prec.P2*pinv(full(dat.A))*(dat.A*(prec.P2\((prec.E2'*mu)/sigma)+xp_bar)-dat.b);
  else
      prev.mu = zeros(n, 1);
      mu      = zeros(n, 1);
   end

% step 3: convergence check
% ATTENTION! x unscaled, p,mu scaled on the dual
  dxp = xp-prev.xp;
  dlambda = lambda-prev.lambda;
  r1 = (1/tau)*(-dxp) + prec.L'*dlambda;
  s1 = (1/sigma)*(-dlambda) + prec.L*dxp;
  if ~isempty(dat.A)
    dmu = mu-prev.mu;
    r2 = dmu;
    s2 = (1/sigma)*(-dmu) + dxp;
    s = [s1; s2];
    r = r1 + r2;
  else
    s = s1;
    r = r1;
  end

  rPrimal = norm(r,1);
  rDual   = norm(s,1);
%   jjj

  % step 4: update stepsize
  tau = theta*tau;
  theta = 1/sqrt(1+tau*2*gamma/lam);
  sigma = sigma / theta;
  
  if RECORD_STATS
     if settings.fun_eval == 1
        stats_plot(iii).f_eval  = obj;
     end
    stats_plot(iii).rDual       = rDual;
    stats_plot(iii).rPrimal     = rPrimal;
    stats_plot(iii).lambda      = lambda;
    stats_plot(iii).x           = prec.Dinv*xp;
    stats_plot(iii).s           = s;
    stats_plot(iii).r           = r;
    stats_plot(iii).tau         = tau;
    stats_plot(iii).sigma       = sigma;
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
   stats_out.lambda  = lambda;
   stats_out.x       = prec.Dinv*xp;
   stats_out.s       = s;
   stats_out.r       = r;
   stats_out.numiter = jjj;
end

sol.problem = 0;
if rDual > settings.dualTol || rPrimal > settings.primalTol
  cprintf('err', 'Failed to reach optimality after %i iterations\n', settings.maxItr);
  sol.problem = 1;
end

sol.x  = prec.Dinv*xp;
sol.lambda  = lambda;
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