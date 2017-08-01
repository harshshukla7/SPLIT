function [sol, stats_out, stats_plot] = AdPrPDA(prob, settings)
%
% Implementation of PDA with optional scaling and varying stepsizes
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
dat.A = sparse(dat.A);
R     = chol(dat.A*dat.A');

% Validate inputs
validate_inputs(prob, settings);

% -----------------------
% Pre-processing
% -----------------------

%% Data
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
else
    settings.tau = 1;
    tau = settings.tau;
end
prec.Dinv = settings.scale_rPrimal; % D^(-1) (preconditioner)
prec.T = settings.scale_rPrimal*settings.scale_rPrimal'; % T = D^(-1)*D^(-T)
if ~isfield(settings, 'scale_rDual') 
    if ~isempty(dat.A)
        settings.scale_rDual = eye(size(L,1)+n); 
    else
        settings.scale_rDual = eye(size(L,1)); 
    end
else
    settings.rho = 1;
    rho = settings.rho;
end
prec.P = settings.scale_rDual'*settings.scale_rDual; % P=E'E
prec.E = settings.scale_rDual; % E (preconditioner)

% % if ~isfield(settings, 'changeBasisPrimal')
% %     settings.changeBasisPrimal = eye(n); 
% % end
% % SP = settings.changeBasisPrimal;
 
if ~isempty(dat.A)
    prec.E1 = prec.E(1:size(L,1),1:size(L,1));
    prec.E2 = prec.E(size(L,1)+1:size(L,1)+n,size(L,1)+1:size(L,1)+n);
    prec.P1 = prec.P(1:size(L,1),1:size(L,1));
    prec.P2 = prec.P(end-n+1:end, end-n+1:end);
    prec.L  = prec.E1*L*prec.Dinv;
    prec.I  = prec.E2*eye(n)*prec.Dinv;
    prec.K   = [prec.L; prec.I];
    Mu       = max(svds(prec.K))^2; % rho(prec.K'*prec.K)^2 = |prec.K'*prec.K|
else 
    prec.E1 = prec.E(1:size(L,1),1:size(L,1));
    prec.E2 = zeros(n,n);
    prec.P1 = prec.P(1:size(L,1),1:size(L,1));
    prec.P2 = zeros(n, n);
    prec.L = prec.E1*L*prec.Dinv;
    prec.K = prec.L;
    Mu     = max(svds(prec.K))^2;
end

%% Stepsizes and adaptation constants
if isfield(settings,'adapt')
    if ~isfield(settings, 'rho') || ~isfield(settings, 'tau')
        if strcmp(settings.adapt, 'Backtrack')
            ranvec = randn(n);
            settings.tau = sqrt(2*norm(ranvec) / norm(prec.K'*prec.K*ranvec));
            tau = settings.tau;
            settings.rho = tau;
            rho = settings.rho;
        elseif strcmp(settings.adapt, 'Lipschitz')
            settings.tau = 0.95/sqrt(Mu);
            tau = settings.tau;
            settings.rho = 0.95/(tau*(Mu));
            rho = settings.rho;
        end
    end
    % adaptive stepsize constants
    settings.alpha = 0.5;
    alpha = settings.alpha;
    settings.Delta = 1.5;
    settings.eta = 0.95;
    settings.s = 1;
    settings.gamma = 0.75;
    gamma = settings.gamma;
    settings.beta = 0.95;
    beta = settings.beta;
elseif ~isfield(settings,'adapt')
        if ~isfield(settings, 'rho') || ~isfield(settings, 'tau')
            settings.tau = 0.95/sqrt(Mu);
            tau = settings.tau;
            settings.rho = 0.95/(tau*(Mu));
            rho = settings.rho;
        end
end
    

%%  Variables initialization and warm-starting
if ~isfield(settings, 'primal_vars'),     settings.primal_vars = zeros(n, 1); end
x = settings.primal_vars; 
x_bar = settings.primal_vars;

if ~isfield(settings, 'dual_vars_p'),     settings.dual_vars_p = zeros(size(L, 1), 1); end
lam = settings.dual_vars_p; 
lamd = prec.E1*lam;

if ~isfield(settings, 'dual_vars_mu'),    settings.dual_vars_mu = zeros(n, 1); end
mu = settings.dual_vars_mu;
mud = prec.E2*mu;

% -----------------------
% Main
% -----------------------
iii = 0;
obj = 0;
mud = x; % left scaled
for jjj = 1 : 2 %settings.maxItr
  iii = iii + 1;
    
  % step 0
  prev.x    = x;
  
  %% step 1: linear system solve
  x = (eye(size(dat.Q))+tau*prec.T*dat.Q) \ (prev.x-tau*prec.T*(L'*lamd+mud+dat.f'))
  x_bar = 2*x - prev.x

  %% step 2: compute projection
  prev.lamd = lamd;
  lamd_tilde = lamd + rho*prec.P1*(L*x_bar + l);
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

  if ~isempty(dat.A) % dynamics update
      prev.mud    = mud;
      mud = rho*prec.P2*dat.A'*(R\(R'\(dat.A*(prec.P2\(mud/rho)+x_bar)-dat.b)));
  else
      prev.mud = zeros(n, 1);
      mud      = zeros(n, 1);
  end
  
% % % %   % step 6: convergence check
% % % %   % ATTENTION! data unscaled, x unscaled, p,mu scaled 
% %   dx = x-prev.x;
% %   dp = -prec.E(1:size(L,1),1:size(L,1))'\(prev.lamd-lamd);
% %   r1 = (1/tau)*(prec.T\(-dx)) + L'*(dp);
% %   s1 = (1/rho)*(prec.P(1:size(L,1),1:size(L,1))\-dp) + L*dx;
% %   if ~isempty(dat.A)
% %     dmu = prec.E(size(L,1)+1:size(L,1)+n,size(L,1)+1:size(L,1)+n)'\(mud-prev.mud);
% %     r2 = dmu;
% %     s2 = (1/rho)*(prec.P(size(L,1)+1:size(L,1)+n,size(L,1)+1:size(L,1)+n)\-dmu) + dx;
% %     s = [s1; s2];
% %     r = r1 + r2;
% %   else
% %     s = s1;
% %     r = r1;
% %   end

%% step 3: convergence check
% ATTENTION! x unscaled, lam, mu scaled on the dual
  dx = x-prev.x;
  dlamd = lamd-prev.lamd;
  r1 = (1/tau)*(prec.Dinv\(-dx)) + prec.Dinv'*(L'*(dlamd));
  s1 = (1/rho)*(prec.E1'\(-dlamd)) + prec.E1*L*dx;
  if ~isempty(dat.A)
    dmud = mud-prev.mud;
    r2 = prec.E2'\dmud;
    s2 = (1/rho)*(prec.E2'\(-dmud)) + prec.Dinv\dx;
    s = [s1; s2];
    r = r1 + r2;
  else
    s = s1;
    r = r1;
  end

  rPrimal = norm(r,1);
  rDual   = norm(s,1);
 jjj
  
  if rDual < settings.dualTol && rPrimal < settings.primalTol
    cprintf([0.1 0.5 0.1], '\n\n>>>>> Stopping on optimality after %i steps <<<<<\n\n', iii);
    break
  end
% % diff = abs(settings.fval-dat.f*x)
% %   if diff <= 1e-2
% %     cprintf([0.1 0.5 0.1], '\n\n>>>>> Stopping on optimality after %i steps <<<<<\n\n', iii);
% %     break
% %   end
  
  %% step 4: backtracking and stepsize adaptation
  if ~isfield(settings,'adapt')
      tau = tau;
      rho = rho;
  elseif strcmp(settings.adapt, 'Backtrack')
      if ~isempty(dat.A)
           b = ( 2*tau*rho * dx' * [L' eye(n)'] * [dlamd; dmud] ) / ...
          ( gamma*rho*dx'*(prec.T\(dx)) + gamma*tau*([dlamd; dmud])'*(prec.P\([dlamd; dmud])) ); % checking condition
      else
           b = ( 2*tau*rho * dx' * L' * dlamd ) / ...
          ( gamma*rho*dx'*(prec.T\(dx)) + gamma*tau*dlamd'*(prec.P\dlamd) ); 
      end
      if b > 1
          tau = beta*tau/b;
          rho = beta*rho/b;
          x = prev.x;
          lamd = prev.lamd;
          mud = prev.mud;
          alpha = settings.alpha;
      elseif (rPrimal > settings.s*rDual*settings.Delta)
          tau = tau / (1-alpha);
          rho = rho * (1-alpha);
          alpha = alpha * settings.eta;
      elseif (rPrimal < settings.s*rDual/settings.Delta)
          tau = tau * (1-alpha);
          rho = rho / (1-alpha);
          alpha = alpha * settings.eta;
      elseif (rPrimal >= settings.s*rDual/settings.Delta) && (rPrimal <= settings.s*rDual*settings.Delta) 
          tau = tau;
          rho = rho;
          alpha = alpha;
      end
  elseif strcmp(settings.adapt, 'Lipschitz')   
      if (rPrimal > settings.s*rDual*settings.Delta)
          tau = tau / (1-alpha);
          rho = rho * (1-alpha);
          alpha = alpha * settings.eta;
      end
  
      if (rPrimal < settings.s*rDual/settings.Delta)
          tau = tau * (1-alpha);
          rho = rho / (1-alpha);
          alpha = alpha * settings.eta;
      end
  
      if (rPrimal >= settings.s*rDual/settings.Delta) && (rPrimal <= settings.s*rDual*settings.Delta) 
          tau = tau;
          rho = rho;
          alpha = alpha;
      end
  end
  
  %% record statistics (optional)
  if RECORD_STATS
     if settings.fun_eval == 1
        stats_plot(iii).f_eval  = obj;
     end
    stats_plot(iii).rDual       = rDual;
    stats_plot(iii).rPrimal     = rPrimal;
    stats_plot(iii).lam         = prec.E1\lamd;
    stats_plot(iii).mu          = prec.E2\mud;
    stats_plot(iii).x           = x;
    stats_plot(iii).s           = s;
    stats_plot(iii).r           = r;
    stats_plot(iii).tau         = tau;
    stats_plot(iii).rho         = rho;
  end 

end

% record outside the loop hence we have only the most recent values
if RECORD_STATS
   stats_out.rDual   = rDual;
   stats_out.rPrimal = rPrimal;
   stats_out.lam     = prec.E1\lamd;
   stats_out.mu      = prec.E2\mud;
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

sol.x  = x;
sol.lam  = prec.E1\lamd;
sol.mu  = prec.E2\mud;
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