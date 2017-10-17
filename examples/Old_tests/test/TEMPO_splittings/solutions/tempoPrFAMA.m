function [sol, stats_out, stats_plot] = tempoPrFAMA(dat, dxs, dus, xinit, settings)
%
% Implementation of AMA with optional scaling, Nesterov acceleration
% and restart
%
% Solves the problem
%   min   x' * Q * x + f' * x + g(L * x - l)
%   s.t.  A * x = b
%
% Structure 'prec' containts the scaling matrix E that should
% always be positive definite diagonal. If not specified set to identity
%
% Variables followed by a 'd' indicate scaling by constraint preconditioner E
%

%% Data
Q = dat.cost.Q;
R = dat.cost.R;
P = dat.cost.P;
A = dat.model.A;
B = dat.model.B;
N = dat.time.N;
n = dat.dim.nx;
m = dat.dim.nu;

%% Scalings
if ~isfield(settings, 'dualPrec') 
    settings.dualPrec = eye(size(dat.L,1)); 
end
prec.P  = settings.dualPrec'*settings.dualPrec; % P=E'E
prec.E  = settings.dualPrec; % E (preconditioner)
L = prec.E*dat.L;   l = prec.E*dat.l;

%% Setup the KKT matrix and pre-solve
%%*******************************************************************%%
%% Build the linear system                                           %%
%% Your QP reads                                                     %%  
%% minimize (z-zs)'*H*(z-zs)                                         %%
%% subject to Cz=d                                                   %%
%%            -Lz+l>=0                                               %%
%% with z=(x,u), zs=(xs,us), H=diag(Q,Q...,P,R,R,,,,R), C:dynamics   %%
%% Formulate the linear system that you will invert in the first step%% 
%% of the algorithm.                                                 %% 
%%*******************************************************************%%

H = 2*blkdiag(kron(eye(N),Q),P,kron(eye(N),R));
tilde.A  = kron(eye(N+1),eye(n)) + [zeros(n,n*(N+1));[kron(eye(N),-A) zeros((N)*n,n)]];
tilde.B  = kron(eye(N),-B); 
C = [tilde.A [zeros(n,m*N);tilde.B]];

M = [ H C'; C 0*speye((N+1)*n,(N+1)*n) ];
b2 = [xinit; zeros(n*N,1)];
b1_const = H*[dxs;dus];
rho = 1; % stepsize

%% Main
beta = 1;
prev.beta = beta;
lam = zeros(size(L,1),1);
z = zeros((N+1)*n+N*m,1);
prev.lam = lam;
iii = 0;
for jjj = 1:settings.maxItr
  iii = iii + 1;
  
  % Step 1: Nesterov acceleration
  beta = (1+sqrt(4*prev.beta^2+1))/2;
  prev.beta = beta;
  hat.lam = lam + (prev.beta-1)*(lam - prev.lam)/beta;

  % Step 2 : solve linear system
  prev.z = z;
  b1_var = b1_const+L'*hat.lam;
  b = [b1_var;b2];
  sol = M \ b;
  z = sol(1:n*(N+1)+m*N);
  
  % Step 3 : compute prox operators
  q =  prec.E\(-L*z + l + (1/rho*hat.lam));
  yd = prec.E*max(q,0);
          
  % Step 4: update the Lagrange multiplier
  prev.lam = lam;
  lam = hat.lam + rho*(-L*z + l - yd);
  
  % Step 5: adaptive restart
  if isfield(settings,'restart')
      % generalized gradient evaluation
      test = (hat.lam-lam)' * (lam-prev.lam); 
      if test > 0 % monotonicity test
          hat.lam = prev.lam; 
          b1_var = b1_const+L'*hat.lam;
          b = [b1_var;b2];
          sol = M \ b;
          z = sol(1:n*(N+1)+m*N);
          q =  prec.E\(-L*z + l + (1/rho*hat.lam));
          yd = prec.E*max(q,0); 
          prev.lam = hat.lam;
          lam = hat.lam + rho*(-L*z + l - yd);    
      end
  end

  % step 4: convergence check
  % L, l scaled, z, lam unscaled, yd scaled
  s = L'*(prev.lam - lam);
  r = -L*z + l - yd;
  
  rPrimal = norm(r,1);
  rDual   = norm(s,1);
  
   
  % record statistics (optional)
 stats_plot(iii).rDual   = rDual;
 stats_plot(iii).rPrimal = rPrimal;
 stats_plot(iii).x       = z;
 stats_plot(iii).y       = prec.E\yd;
 stats_plot(iii).lam     = lam;
 stats_plot(iii).s       = s;
 stats_plot(iii).r       = r;
 stats_plot(iii).rho     = rho;
  
  if rDual < settings.dualTol && rPrimal < settings.primalTol
    cprintf([0.1 0.5 0.1], '\n\n>>>>> Stopping on optimality after %i steps <<<<<\n\n', iii);
    break
  end

end
  
% record outside the loop hence we have only the most recent values
   stats_out.rDual   = rDual;
   stats_out.rPrimal = rPrimal;
   stats_out.x       = z;
   stats_out.y       = prec.E\yd;
   stats_out.lam     = lam;
   stats_out.s       = s;
   stats_out.r       = r;
   stats_out.rho     = rho;
   stats_out.numiter = jjj;

if rDual > settings.dualTol || rPrimal > settings.primalTol
  cprintf('err', 'Failed to reach optimality after %i iterations\n', settings.maxItr);
  sol.problem = 1;
end

sol.lam = lam;
sol.x = z;
sol.y = prec.E\yd;

end