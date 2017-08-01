function [sol, stats_out, stats_plot] = AdPrADMM(prob, settings)
%
% Quick implementation of ADMM with optional scaling and adaptive stepsize
%
% Solves the problem
%   min 0.5 x' * Q * x + f' * x + c+ sum_i Fi(Li * x + li)
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
if ~isfield(settings, 'rho'),                  settings.rho = 1; end
if ~isfield(settings, 'relaxation'),           settings.relaxation = 1; end
if ~isfield(settings, 'dualTol'),              settings.dualTol = 1e-4; end
if ~isfield(settings, 'primalTol'),            settings.primalTol = 1e-4; end
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

% Relaxation constant for adaptive stepsize
rho   = settings.rho;
if isfield(settings, 'adapt')
    alpha = 1;
    gamma = (1+sqrt(5)) / 2 - eps;
else
    alpha = settings.relaxation;
    gamma = 1;
end

% -----------------------
% Pre-processing
% -----------------------
%% data
n = size(dat.Q,1);
m = size(dat.A,1);


%%


% Form the augmented lagrangian
% min 0.5*x'*Q*x + f'*x + sum_i Fi(y_i) + rho/2*sum_i ||Li*x + li - y_i + lam_i||
%      A*x = b

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
prec.Dinv = sparse(settings.scale_rPrimal);
prec.T = settings.scale_rPrimal'*settings.scale_rPrimal;
if ~isfield(settings, 'scale_rDual') 
    settings.scale_rDual = eye(size(L,1)); 
end
prec.P  = sparse(settings.scale_rDual'*settings.scale_rDual); % P=E'E
prec.E  = sparse(settings.scale_rDual); % E (preconditioner)

% % if ~isfield(settings, 'changeBasisPrimal')
% %     settings.changeBasisPrimal = eye(n); 
% % end
% % SP = settings.changeBasisPrimal;

%% Setup the KKT matrix and pre-solve
Q = dat.Q;
LL = zeros(n,n);
sum_l = zeros(n,1);
sum_L = [];
L = []; l = []; 
nProxVars = zeros(nProx,1);
kk = 1;
for i = 1:nProx
  proxInd{i}   = kk + [0:size(prox(i).L,1)-1];
  ind          = proxInd{i};
  sum_L        = [sum_L prox(i).L'*prec.E(ind,ind)'];
  sum_l        = sum_l + rho*prox(i).L'*(prec.E(ind,ind)'*prec.E(ind,ind))*prox(i).l;
  l            = [l;prec.E(ind,ind)*prox(i).l]; % ATTENTION! L and l are left scaled
  L            = [L;prec.E(ind,ind)*prox(i).L];
  Q            = Q + rho*prox(i).L'*(prec.E(ind,ind)'*prec.E(ind,ind))*prox(i).L;
  LL           = LL + (prec.E(ind,ind)*prox(i).L)'*(prec.E(ind,ind)*prox(i).L);
  kk           = kk + size(prox(i).L,1);
  nProxVars(i) = size(prox(i).L,1);
end

%%

LL = full(LL);
dat.A = full(dat.A);
proxWeight = [prob.prox.weight]; % Weights for the prox functions

%% Simultaneous diagonalization of KKT system
if isfield(settings, 'adapt')
    if (min(eig(LL))>0) && (isempty(dat.A))
        H = dat.Q;   M = LL;
        Lchol = chol(M,'lower');
        C = Lchol \ (H*inv(Lchol'));
        [V,D1] = qdwheig(C);
        cas = 1;
    elseif (min(eig(dat.Q))>0) && (isempty(dat.A))
        H = LL;   M = dat.Q;
        Lchol = chol(M,'lower');
        C = Lchol \ (H*inv(Lchol'));
        [V,D1] = qdwheig(C);
        cas = 2;
    else
        error('Conditions for simultaneous diagonalization not satisfied. Adaptive ADMM not adviced - increased cost per iteration.')
        settings.adapt = 'no';
    end
    X = Lchol' \ V;  
    % Permute D - was not working otherwise
    D1 = blkdiag(diag(nonzeros(diag(X'*H*X))),zeros(size(D1,1)-size(nonzeros(diag(X'*H*X)),1),size(D1,1)-size(nonzeros(diag(X'*H*X)),1)));   
    D1 = sparse(D1);
else
    % Solution to KKT system
    KKT = [Q dat.A'; dat.A zeros(m)];
    [Lkkt,Dkkt,Pkkt] = ldl(KKT);
    
    rhsKKTconst = [-dat.f' - sum_l; dat.b];
    rhsKKTvar   = [-rho*sum_L;sparse(length(dat.b),size(sum_L,2))];
    const_vector = Pkkt*(Lkkt'\(Dkkt\(Lkkt\(Pkkt'* rhsKKTconst))));
    matrix = Pkkt*(Lkkt'\(Dkkt\(Lkkt\(Pkkt' * rhsKKTvar))));
end

%%  Variables initialization and warm-starting
if ~isfield(settings, 'primal_vars_x'),     settings.primal_vars_x = zeros(n, 1); end
x = settings.primal_vars_x; 
if ~isfield(settings, 'primal_vars_y'),    settings.primal_vars_y = zeros(size(L,1),1); end
y  = settings.primal_vars_y;
yd = prec.E*y;
if ~isfield(settings, 'dual_vars_lam'),     settings.dual_vars_lam = zeros(size(L,1),1); end
lam = settings.dual_vars_lam;         

% -----------------------
% Main
% -----------------------
iii = 0;
for jjj = 1 : settings.maxItr
  iii = iii + 1;

  %% step 1: linear system solve
  if isfield(settings, 'adapt')
      sum_l = zeros(n,1);
      for i = 1:nProx
        sum_l     = sum_l + rho*prox(i).L'*(prec.E(ind,ind)'*prec.E(ind,ind))*prox(i).l;
      end
      rhsKKTconst = [-dat.f' - sum_l; dat.b];
      rhsKKTvar   = [-rho*sum_L; sparse(length(dat.b),size(sum_L,2))];

      if cas==1
          rho1=1; rho2=rho;
      elseif cas==2
          rho1=rho; rho2=1;
      end
      x = (X * ((rho1*D1+rho2*eye(n)) \ X'))* (rhsKKTconst + rhsKKTvar*(lam/rho - yd));
  else
      %disp('sum_l is')
      %sum_l
      t = const_vector + matrix * (lam/rho - yd); 
      x = t(1:n);
      %rho
  end
  
  %% step 2 : compute prox operators
  
  Lx_relax = alpha*L*x - (1-alpha)*(-yd+l);
  q =  prec.E \ (Lx_relax + l + 1/rho * lam);
  prev.yd = yd;
  for i = 1:nProx
    ind = proxInd{i};
    yd(ind) = diag(prec.E(ind,ind)).*prox(i).func(q(ind), 1/rho*proxWeight(i)); 
  end
          
  %% step 3: Update the lagrange multiplier
  lam = lam + gamma*rho*(Lx_relax + l - yd);
  
  
% %   % function evaluation
% %   if settings.fun_eval == 1
% %       obj = .5*(SP\x)'*dat.Q*(SP\x) + dat.f*(SP\x);
% %   end
  
% %   %% Convergence check on original
% %   % L, l scaled, z, lam unscaled, yd scaled
% %   s = rho*(prec.E\L)'*(prec.E\(prev.yd-yd));
% %   r = prec.E\(L*x + l - yd);

  %% step 4: convergence check
  % L, l scaled, z, lam unscaled, yd scaled 
  s = rho*(L*prec.Dinv)'*(prev.yd-yd);
  r = L*x + l - yd;
  
  rPrimal = norm(r,1);
  rDual   = norm(s,1);
  
  
  %% step 5: adapt stepsize
  if isfield(settings, 'adapt')
      if rPrimal > 2*rDual; rho = 2*rho;
      elseif rDual > 2*rPrimal; rho = 0.5*rho;
      end
  end
  
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
%    stats_out.KKT     = KKT;
end

sol.problem = 0;
if rDual > settings.dualTol || rPrimal > settings.primalTol
  cprintf('err', 'Failed to reach optimality after %i iterations\n', settings.maxItr);
  sol.problem = 1;
end

sol.lam = lam;
sol.x = x;
sol.y = prec.E\yd;

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