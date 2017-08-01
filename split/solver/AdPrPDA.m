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
sumMu = 0;
for i = 1:nProx
  sumMu        = sumMu + max(svds(prox(i).L))^2; 
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
    if (~isempty(prob.dat.Q))&&(nnz(prob.dat.Q))  % existence of smooth term changes stepsize
        settings.tau = 1/max(eig(dat.Q));
        tau = settings.tau;
        settings.rho = (1-tau*max(eig(dat.Q))/2)^2 / (tau*Mu);
        rho = settings.rho;
        % over-relaxation
        kappa = (1/max(eig(dat.Q)))*(1/tau-rho*Mu);
        delta = max(4*kappa/(2*kappa+1),min(3/2,1/2+kappa));
        theta = delta;
    else
        settings.tau = 1/sqrt(Mu);
        tau = settings.tau;
        settings.rho = 1/(tau*Mu);
        rho = settings.rho;
        % over-relaxation
        theta = 1.8;
    end
    if isfield(settings,'adapt') % No prec, Yes adapt
        if any(dat.Q)
            settings.tau = 0.05/max(eig(dat.Q));
            tau = settings.tau;
            settings.rho = (1-tau*max(eig(dat.Q))/2)^2 / (tau*sumMu);
            rho = settings.rho;
        else
            settings.tau = 1/sqrt(sumMu);
            tau = settings.tau;
            settings.rho = tau;
            rho = settings.rho;
        end
        theta = 1;
        % adaptive stepsize constants
        settings.alpha = 0.5;
        alpha = settings.alpha;
        settings.Delta = 1.5;
        settings.eta = 0.95;
        settings.s = 1;
    end
end
prec.Dinv = sparse(settings.scale_rPrimal); % D^(-1) (preconditioner)
prec.D = sparse(inv(settings.scale_rPrimal)); % D
prec.T = sparse(settings.scale_rPrimal*settings.scale_rPrimal'); % T = D^(-1)*D^(-T)
prec.P  = settings.scale_rDual'*settings.scale_rDual; % P=E'E
prec.E  = sparse(settings.scale_rDual); % E (preconditioner)
prec.L  = prec.E*L*prec.Dinv;

if isfield(settings, 'precondition') 
    Mu    = settings.Mu;
    theta = 1;
    %% Stepsizes and adaptation constants
    if isfield(settings,'adapt') % Yes prec, Yes adapt
        if any(dat.Q)
            settings.tau = 1/(max(svds(prec.T))*max(eig(dat.Q)));
            tau = settings.tau;
            settings.rho = (1-tau*max(svds(prec.T))*max(eig(dat.Q))/2)^2 / (tau*Mu);
            rho = settings.rho;
        else
            settings.tau = 1/sqrt(Mu);
            tau = settings.tau;
            settings.rho = tau;
            rho = settings.rho;
        end
        % adaptive stepsize constants
        settings.alpha = 0.5;
        alpha = settings.alpha;
        settings.Delta = 2;
        settings.eta = 0.95;
        settings.s = 1;
    elseif ~isfield(settings,'adapt') % Yes prec, No adapt
        rho = 1;
        tau = 1;
    end
end
    

%%  Variables initialization and warm-starting
if ~isfield(settings, 'primal_vars_x'),     settings.primal_vars_x = zeros(n, 1); end
x = settings.primal_vars_x; 
% x_bar = settings.primal_vars;

if ~isfield(settings, 'dual_vars_lam'),     settings.dual_vars_lam = zeros(size(L, 1), 1); end
lam = settings.dual_vars_lam; 
lamd = prec.E*lam;

% -----------------------
% Main
% -----------------------
iii = 0;
obj = 0;
for jjj = 1 : settings.maxItr
  iii = iii + 1;
  
  %% step 1: linear system solve
  
  
  prev.x = x;
  y = (eye(size(dat.Q))-tau*prec.T*dat.Q)*x-tau*prec.T*(L'*lamd+dat.f');
  if ~isempty(dat.A) % dynamics update
      p = y + dat.A'*((dat.A*dat.A')\(dat.b-dat.A*y));
  else
      p = y;
  end
  
  x_bar = 2*p-prev.x;
   
  x = x + theta*(p-x);

  
  %% step 2: compute projection
  prev.lamd = lamd;
  lamd_tilde = lamd + rho*prec.P*(L*x_bar + l);
  
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
      mud(ind) = prox(i).conjFunc(lamd_tilde(ind), proxWeight(i));
      prev.lamd(ind) = lamd(ind);
      lamd(ind) = lamd(ind) + theta*(mud(ind)'-lamd(ind));
  end
  
  
  
  
% %   % step 6: convergence check
% %   % ATTENTION! data unscaled, x unscaled, p,mu scaled 
  dx = x-prev.x;
  dlamd = lamd-prev.lamd;
  r = (dat.Q-(1/tau)*prec.D)*dx + prec.Dinv'*(L'*(dlamd));
  s = (1/rho)*(prec.E'\(-dlamd)) + prec.E*L*dx;

% % 
% % %% step 3: convergence check
% % % ATTENTION! x unscaled, lam, mu scaled on the dual
% %   dx = x-prev.x;
% %   dlamd = lamd-prev.lamd;
% %   r1 = (1/tau)*(prec.Dinv\(-dx)) + prec.Dinv'*(L'*(dlamd));
% %   s1 = (1/rho)*(prec.E1'\(-dlamd)) + prec.E1*L*dx;
% %   if ~isempty(dat.A)
% %     dmud = mud-prev.mud;
% %     r2 = prec.E2'\dmud;
% %     s2 = (1/rho)*(prec.E2'\(-dmud)) + prec.Dinv\dx;
% %     s = [s1; s2];
% %     r = r1 + r2;
% %   else
% %     s = s1;
% %     r = r1;
% %   end
% % 
  rPrimal = norm(r,1);
  rDual   = norm(s,1);
  %jjj
  
  %% step 4: backtracking and stepsize adaptation
  if ~isfield(settings,'adapt')
      tau = tau;
      rho = rho;
  elseif strcmp(settings.adapt, 'Lipschitz')   
      if (rPrimal > settings.s*rDual*settings.Delta)
          prev.tau = tau;
          prev.rho = rho;
          prev.alpha = alpha;
          tau = tau / (1-alpha);
          rho = rho * (1-alpha);
          alpha = alpha * settings.eta;
          tau = tau * (1-alpha);
          rho = rho / (1-alpha);
          alpha = alpha * settings.eta;
          % check convergence condition
          if ( sqrt(tau*rho) >= (1-max(eig(dat.Q))/2*max(svds(prec.T))*tau+eps)/sqrt(Mu) )
              tau = prev.tau;
              rho = prev.rho;
              alpha = prev.alpha;
          end
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

sol.x  = x;
sol.y  = y;
sol.lam  = prec.E\lamd;
sol.lam_noscale = lamd;
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