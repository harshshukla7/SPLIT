% Everything before the IO line is executed on object creation
% This means that none of the parameters, constants, etc have been defined yet
% Use this block to define the inputs and outputs

% Solves problems of the form
%         min z'*z + 2*h'*z;
%         K*z <= k;

assert(isempty(cdr.A), 'VM FAMA can currently only handle problems with no equality constraints')
assert(length(cdr.prob.prox) == 1 && strcmp(cdr.prob.prox(1).typeName, 'nonNegative'), ...
  'VM FAMA can currently only handle quadratic programs')
assert(norm(2*eye(size(cdr.Q,1)) - cdr.Q) < 1e-6, 'Q must be the identity')

% Parameters available at compile-time and run-time
cdr.addConst('rho',      'Penalty parameter');
cdr.addConst('itr_conv', 'Number of iterations between convergence checks')
cdr.addConst('VM',       'Variable metric weighting')

cdr.addConst('SP',       'Primal residual scaling')
% cdr.addConst('SD',       'Dual residual scaling')

% Arguments to the generated function
% cdr.addInput('iterationLimit', [1 1], 'Run-time upper limit on number of iterations')
% cdr.addInput('runTimeDualTol', [1 1], 'Run-time required tolerance on dual')

% Outputs of the generated function
cdr.addOutput('y_',      'Copy of the primal variable');
cdr.addOutput('lambda_', 'Dual variable');
cdr.addOutput('itr',     'Number of iterations');
cdr.addOutput('tm',      'Computation time');
cdr.addOutput('rDual',   'Dual residual');

% Parameters available at compile-time, but not run-time
cdr.addParameter('KKTSolveMode',       '''LU'',''Invert KKT''')
cdr.addParameter('MatVec_SparseLimit', 'Percent sparsity below which matrix-vector products use code-gen')

% >>>>>>>>>>> IO >>>>>>>>>>>
%
% VM FAMA template
%
% Solves problems of the form
%         min z'*z + 2*h'*z;
%         K*z <= k;
%
% >>>>>>>>>>> HEADER >>>>>>>>>>>

% >>>>>>>>>>> RUN >>>>>>>>>>>

% Allocate memory for internal variables
% Note that output variables must have different names from internal
% variables. Use x_ as output, and x as internal - then copy x_ to x at the
% end of the algorithm
cdr.addVar('x',       size(cdr.L,2),1); % Size argument can be a string, evaluated at compile-time
cdr.addVar('y',       length(cdr.l),1); % Size argument can be a constant
cdr.addVar('lambda',  length(cdr.l),1);
cdr.addVar('lam_hat', length(cdr.l),1);

% Compute the parametric solution
cdr.genComputeParametric('par');

%% >>>>>>>>>>> COPY >>>>>>>>>>>
rDual   = 0;
rPrimal = 0;
itr     = 0;
beta_k  = 1;

for iii = 1:maxItr
  % Stop on run-time iteration limit
  if iii >= iterationLimit
    break
  end
  
  prev_lam = lambda;
  prev_y      = y;
  
  % Step 1 : solve linear system
  % >>>>>>>>>>> RUN >>>>>>>>>>
  
  % workLp = L' * lam_hat
  cdr.genProd('workLp', -cdr.L', '(-L'')', 'lam_hat');
  
  % Pre-solve the KKT system
  % x = Q \ (-f + workLp)
  %   cdr.mldivide(cdr.Q, '(-f + workLp)', 'x', cdr.KKTSolveMode);
  
  % >>>>>>>>>>> COPY >>>>>>>>>>
  x = (-f + workLp) ./ diag(Q); % Q is diagonal and invertable
  
  % Update the lagrange multiplier
  prev_lam = lambda;
  lambda      = min(0, lam_hat + VM.*(L*x + l));
  
  % Nesterov acceleration:
  beta = beta_k;
  beta_k = (1+sqrt(4*beta*beta+1))/2;
  lam_hat = lambda + (beta-1)*(lambda - prev_lam)/beta_k;
  
  % Convergence checks
  if mod(iii, itr_conv) == 0
    % Convergence checks
    s = -SP.*(L'*(prev_lam - lambda));
    rDual   = s'*s;
    if rDual < (dualTol*dualTol) && rDual < (runTimeDualTol*runTimeDualTol)
      break
    end
  end
end

% Set the output
itr     = iii;
lambda_ = lambda;
x_      = x;
y_      = y;
rDual   = sqrt(rDual);

%% >>>>>>>>>>> RUN >>>>>>>>>>>
% Must call this to terminate the function (and compute the computation time)
% All addFunction calls must come before this point
cdr.genFunctionEnd();

