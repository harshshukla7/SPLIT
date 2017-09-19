function [sol, stats] = chambolle_prec(prob, settings)
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
%    prox.type :
%      - 'ind_nonnegative' : indicator for y >= 0
%      - 'ind_norm_i_cone' : indicator for norm(y(1:end-1), i) <= y(end)
%      - 'Li'              : norm(y, i)
%      - 'ind_quad'        : indicator for y'*y <= 1
%    prox.L, l Data matrices
%
% settings.rho = stepsize
%

if nargin < 2
  settings = [];    
end

% Convergence statisctics requested
RECORD_STATS = false;
if nargout > 1 
  RECORD_STATS = true;
end

dat   = prob.dat;
prox  = prob.prox;
conj  = prob.conj;
nProx = length(prox);

Acal = [];
for i = 1:nProx
    Acal = [Acal; prox(i).L];
end
Acal = [Acal; dat.A];
% radius = sqrt(max(eigs(Acal'*Acal)));

% % isdiag = @(A) nnz(A)==nnz(diag(A));

% Validate inputs
validate_inputs(prob, settings);

% Diagonal preconditioning and relaxation
alpha = 1;
sumsigma = zeros(size(Acal,1),1);
sumtau = zeros(size(Acal,2),1);
for j = 1:size(Acal,2)
 for i = 1:size(Acal,1)
     sumtau(j) = sumtau(j) + abs(Acal(i,j))^(2-alpha);
     sumsigma(i) = sumsigma(i) + abs(Acal(i,j))^alpha;
 end
end

Tau = diag(1./sumtau);
Sigma = diag(1./sumsigma);
theta = 1;

if ~isfield(settings, 'terminationCheckFreq'), settings.terminationCheckFreq = 10; end
if ~isfield(settings, 'Sigma'),                settings.Sigma = diag(Sigma); end
if ~isfield(settings, 'Tau'),                  settings.Tau = diag(Tau); end
if ~isfield(settings, 'relaxation'),           settings.relaxation = theta; end
if ~isfield(settings, 'pre-condition'),        settings.alpha = alpha; end
if ~isfield(settings, 'dualTol'),              settings.dualTol = 1e-3; end
if ~isfield(settings, 'primalTol'),            settings.primalTol = 1e-3; end
if ~isfield(settings, 'maxItr'),               settings.maxItr = 2e3; end


% -----------------------
% Pre-processing
% -----------------------

n = size(dat.Q,1);
m = size(dat.A,1);

% Setup the KKT matrix and pre-solve
Q = dat.Q;
% sum_l = zeros(n,1);
sum_L = [];
L = []; l = []; 
nProxVars = zeros(nProx,1);
kk = 1;
for i = 1:nProx
%   sum_l        = sum_l + rho*prox(i).L'*prox(i).l;
  sum_L        = [sum_L prox(i).L'];
  l            = [l;prox(i).l];
  L            = [L;prox(i).L];
  proxInd{i}   = kk + [0:size(prox(i).L,1)-1];
  kk           = kk + size(prox(i).L,1);
  nProxVars(i) = size(prox(i).L,1);
end


% Solution to Linear System:
%   [x;dual] = iKKT * (rhsKKTconst + rhsKKTvar * q)
% where q = lambda - y
LS         = Q + Tau\eye(n);
rhsLSconst = -dat.f';
rhsLSvar   = [-sum_L Tau\eye(n) -dat.A'];

% Primal-dual formulation
% min max 0.5*x'*Q*x + f'*x + <A*x-b,nu> + sum_i ( <Li*x+li,pi> - Fi*(pi) )
%      

% Matrices to store the solution
z      = zeros(n,1); 
bar.z  = zeros(n,1);
p      = zeros(sum(nProxVars),1);
nu     = zeros(m,1);
q      = zeros(size(L,1),1);
    

iii = 0;
for jjj = 1:settings.maxItr
  iii = iii + 1;
%   jjj
  % Step 1 : compute prox operators
  % q = p + sigma(Lbar.z+l)
  old.p  = p;
  old.nu = nu;
  q = p + L*bar.z+l;
  for i = 1:nProx
    ind = proxInd{i}; % identifies the length of pi
    conj(i).scaling = Sigma(ind,ind)*conj(i).scaling;  % XXX careful at prox_norm!!!
    p(ind) = conj(i).conjFunc(q(ind)); 
  end
  nu = nu + Sigma(ind(end)+1:end,ind(end)+1:end)*(dat.A*bar.z-dat.b);
  
  % Step 2 : solve linear system
  old.z = z;
  z = LS \ (rhsLSconst + rhsLSvar*[p;z;nu]);

  % Compute over-relaxation
  bar.z = z + theta*(z-old.z);
%   Lx_hat = L*x;
  
  % Convergence checks
  s  = z-old.z;
  r  = Acal*(z-bar.z) + Sigma \ [p-old.p; nu-old.nu];
  
  rDual   = norm(s);
  rPrimal = norm(r);
    
  if rDual < settings.dualTol && rPrimal < settings.primalTol
    cprintf([0.1 0.5 0.1], '\n\n>>>>> Prec. Chambolle: Stopping on optimality after %i steps <<<<<\n\n', iii);
    break
  end
  
  if RECORD_STATS
    stats(iii).rDual    = rDual;
    stats(iii).rPrimal  = rPrimal;
    stats(iii).z        = z;
    stats(iii).p        = p;
    stats(iii).nu       = nu;
    stats(iii).r        = r;
    stats(iii).s        = s;
    stats(iii).Tau      = diag(Tau);
    stats(iii).Sigma    = diag(Sigma);
  end
end

sol.problem = 0;
if rDual > settings.dualTol || rPrimal > settings.primalTol
  cprintf('err', 'Failed to reach optimality after %i iterations\n', settings.maxItr);
  sol.problem = 1;
end

sol.z  = z;
sol.p  = p;
sol.nu = nu;

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
