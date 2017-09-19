function [sol, stats_out, stats_plot] = PDA(prob, settings)
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

% Extract Lx + l in matrix form
L = []; l = []; 
nProxVars = zeros(nProx,1);
kk = 1;
Mu1 = 0;
for i = 1:nProx
  l            = [l;prox(i).l];
  L            = [L;prox(i).L];
  Mu1          = Mu1 + max(svds(L))^2;
  proxInd{i}   = kk + [0:size(prox(i).L,1)-1];
  kk           = kk + size(prox(i).L,1);
  nProxVars(i) = size(prox(i).L,1);
end
proxWeight = [prob.prox.weight];

% Bound on K operator 
if ~isempty(dat.A)
    Mu2 = 1;
    K   = [L; eye(n)];
else 
    Mu2 = 0;
    K = [L; zeros(n,n)];
end 

% Stepsizes
if ~isfield(settings, 'tau') && ~isfield(settings, 'sigma') 
    tau = 1e-2;
    sigma = 0.5 * (2-tau)^2 / (4*tau*(Mu1+Mu2));
end

% Scalings
if ~isfield(settings, 'scale_rPrimal')
    settings.scale_rPrimal = eye(n); 
else
    tau = 1;
end
prec.Dinv = settings.scale_rPrimal;
prec.T = settings.scale_rPrimal'*settings.scale_rPrimal;
if ~isfield(settings, 'scale_rDual') 
    if ~isempty(dat.A)
        settings.scale_rDual = eye(size(L,1)+n); 
    else
        settings.scale_rDual = eye(size(L,1)); 
    end
else
    sigma = 1;
end
prec.P = settings.scale_rDual'*settings.scale_rDual; % P=E'E
prec.E = settings.scale_rDual; % E (preconditioner)

if ~isfield(settings, 'changeBasisPrimal')
    settings.changeBasisPrimal = eye(n); 
end
SP = settings.changeBasisPrimal;

%  Variables initialization and warm-starting + export
if ~isfield(settings, 'primal_vars'),     settings.primal_vars = zeros(n, 1); end
x = settings.primal_vars; 
x_bar = settings.primal_vars;

if ~isfield(settings, 'dual_vars_p'),     settings.dual_vars_p = zeros(size(L, 1), 1); end
p = settings.dual_vars_p;
pd = p; % left scaled

if ~isfield(settings, 'dual_vars_mu'),    settings.dual_vars_mu = zeros(size(dat.A, 1), 1); end
mu = settings.dual_vars_mu;

% Main
iii = 0;
obj = 0;
mud = x; % left scaled
for jjj = 1:settings.maxItr
  iii = iii + 1;
    
  % step 0
  x_prev    = x;
  
  % step 1
  x = (eye(size(dat.Q))+tau*prec.T*dat.Q) \ (x_prev-tau*prec.T*(K'*[pd;mud]+dat.f'));
  x_bar = 2*x - x_prev;
     
  pd_tilde = pd + sigma*prec.P(1:size(L,1),1:size(L,1))*(L*x_bar + l);
%     nu_prev   = nu;
%     nu = nu + sigma*(dat.A*x_bar - dat.b);
  % step 2: compute projection
  pd_prev    = pd;
%     p_nu_prev = [p_prev; nu_prev];
  for i = 1:nProx
    % special cases: normBall, ellipse - I need to scale the proxWeight vector
    % Perhaps it is not the best solution, but then I need to change the 
    % interface, i.e. add one more field to prob structure e.g. sigma,
    % however, it is a specific case occuring only in chambolle
    if strcmp(char(prox(i).typeConj), 'normProx')
%          proxWeight(i) = sigma*prec.P( proxInd{i}, proxInd{i}) * prox(i).dat.c;
         proxWeight(i) = sigma * prox(i).dat.c;
    end
    if strcmp(char(prox(i).typeConj),  'ellipseConj')
%          proxWeight(i) = sigma*prec.P( proxInd{i}, proxInd{i}) * sqrt(prox(i).dat.c); 
         proxWeight(i) = sigma * sqrt(prox(i).dat.c); 
    end
      ind = proxInd{i};
      pd(ind) = prox(i).conjFunc(pd_tilde(ind), proxWeight(i)); 
  end
  
  if ~isempty(dat.A)
      mud_prev    = mud;
      mud = prec.P(end-n+1:end, end-n+1:end)*pinv(full(dat.A))*(dat.A*(prec.P(end-n+1:end, end-n+1:end)\mud+x_bar)-dat.b);
  else
      mud_prev = zeros(n, 1);
      mud      = zeros(n, 1);
  end
  
  % step 6: convergence check
  % ATTENTION! data unscaled, x unscaled, p,mu scaled 
  dx = x-x_prev;
  r1 = (1/tau)*prec.T\((-dx)) + L'*(prec.E(1:size(L,1),1:size(L,1))'\(pd - pd_prev));
  s1 = (1/sigma)*prec.P(1:size(L,1),1:size(L,1))\((prec.E(1:size(L,1),1:size(L,1))'\(pd_prev-pd))) + L*(dx);
  if ~isempty(dat.A)
    r2 = prec.E(size(L,1)+1:size(L,1)+n,size(L,1)+1:size(L,1)+n)'\(mud-mud_prev);
    s2 = (1/sigma)*prec.P(size(L,1)+1:size(L,1)+n,size(L,1)+1:size(L,1)+n)\(...
    (prec.E(size(L,1)+1:size(L,1)+n,size(L,1)+1:size(L,1)+n)\(mud_prev-mud))) + dx;
    s = [s1; s2];
    r = r1 + r2;
  else
    s = s1;
    r = r1;
  end


  rPrimal = norm(r);
  rDual   = norm(s);
%   jjj
  
  if RECORD_STATS
     if settings.fun_eval == 1
        stats_plot(iii).f_eval  = obj;
     end
    stats_plot(iii).rDual       = rDual;
    stats_plot(iii).rPrimal     = rPrimal;
    stats_plot(iii).p           = p;
    stats_plot(iii).x           = x;
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
   stats_out.p       = p;
   stats_out.x       = x;
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