function [sol, stats_out, stats_plot, data_to_export] = testPDA(prob, settings)
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
if ~isempty(dat.A)
    K = [L; eye(n)];
else
    K = L;
end

% Bound on K operator
Mu = max(svds(K));  

% stepsizes
if ~isfield(settings, 'tau')           
    settings.tau = 1/sqrt(Mu);
    Tau = settings.tau*eye(size(K,2));
else
    Tau = settings.tau;
end
if ~isfield(settings, 'sigma')         
    settings.sigma = settings.tau;
    Sigma = settings.sigma*eye(size(K,1)); % tau and sigma set according to Theorem 2 of paper
else
    Sigma = settings.sigma;
end

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

if ~isfield(settings, 'dual_vars_p'),     settings.dual_vars_p = zeros(size(L, 1), 1); end
p = settings.dual_vars_p;    
data_to_export.p = p;

if ~isfield(settings, 'dual_vars_nu'),    settings.dual_vars_nu = zeros(size(dat.A, 1), 1); end
nu = settings.dual_vars_nu;
data_to_export.nu = nu;

% Scaling of the residuals
if ~isfield(settings, 'scale_rPrimal'),   settings.scale_rPrimal = eye(n); end
SP = settings.scale_rPrimal;
if ~isfield(settings, 'scale_rDual'),    settings.scale_rDual = eye(size(L,1)); end
SD = settings.scale_rDual;
invSP = inv(SP);
invSD = inv(SD);

% Main
iii = 0;
obj = 0;
y = x;
for jjj = 1:settings.maxItr
  iii = iii + 1;
    
  % step 0
  x_prev    = x;
  
  % step 1
  x = (eye(size(dat.Q))+Tau*dat.Q) \ (x_prev-Tau*(L'*p+dat.f'));
  x_bar = 2*x - x_prev;
     
  p_tilde = p + Sigma(1:size(L,1),1:size(L,1))*(L*x_bar + l);
%     nu_prev   = nu;
%     nu = nu + sigma*(dat.A*x_bar - dat.b);
  % step 2: compute projection
  p_prev    = p;
%     p_nu_prev = [p_prev; nu_prev];
  for i = 1:nProx
    % special cases: normBall, ellipse - I need to scale the proxWeight vector
    % Perhaps it is not the best solution, but then I need to change the 
    % interface, i.e. add one more field to prob structure e.g. sigma,
    % however, it is a specific case occuring only in chambolle
    if strcmp(char(prox(i).typeConj), 'normProx')
         proxWeight(i) = Sigma( proxInd{i}, proxInd{i}) * prox(i).dat.c;
    end
    if strcmp(char(prox(i).typeConj),  'ellipseConj')
         proxWeight(i) = Sigma( proxInd{i}, proxInd{i}) * sqrt(prox(i).dat.c); 
    end
      ind = proxInd{i};
      p(ind) = prox(i).conjFunc(p_tilde(ind), proxWeight(i)); 
  end
  
% %   % function evaluation
% %   if settings.fun_eval == 1
% %       obj = .5*(invSP*x)'*dat.Q*(invSP*x) + dat.f*(invSP*x);
% %   end
  
  % step 6: convergence check
%   scaledK = theta*-K;  
  r  = Tau\(invSP*(x_prev-x)) - (L*SP)'*(p_prev - p); 
  s  = Sigma\SD'*([p_prev]-[p]) + invSD*K*(x - x_prev);

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
    stats_plot(iii).s           = s;
    stats_plot(iii).r           = r;
    stats_plot(iii).Tau         = Tau;
    stats_plot(iii).Sigma       = Sigma;
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