% Everything before the IO line is executed on object creation
% This means that none of the parameters, constants, etc have been defined yet
% Use this block to define the inputs and outputs

% Parameters available at compile-time and run-time
cdr.addConst('rho',      'Penalty parameter');
cdr.addConst('itr_conv', 'Number of iterations between convergence checks')

% Arguments to the generated function
% cdr.addInput

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
% ADMM template
%
% Add any custom comments here
%
% This whole block must be comments (no spaces), so that the automatically
% generated comments appear when 'help xxx' is typed at the command line
%
% >>>>>>>>>>> HEADER >>>>>>>>>>>

% >>>>>>>>>>> RUN >>>>>>>>>>>

% Allocate memory for internal variables
% Note that output variables must have different names from internal
% variables. Use x_ as output, and x as internal - then copy x_ to x at the
% end of the algorithm
cdr.addVar('x',     'size(L,2)',1);    % Size argument can be a string, evaluated at compile-time
cdr.addVar('y',     length(cdr.l),1);  % Size argument can be a constant
cdr.addVar('lambda',length(cdr.l),1);

% Compute the parametric solution
cdr.genComputeParametric('par');

%% >>>>>>>>>>> COPY >>>>>>>>>>>
rDual   = 0;
rPrimal = 0;
itr     = 0;

for iii = 1:maxItr
  prev_lambda = lambda;
  prev_y      = y;
  
  % Solve x+
  % >>>>>>>>>>> RUN >>>>>>>>>>
  
  % workLp = L' * ( rho * (-l+y-lambda) )
  cdr.genProd('workLp', cdr.L', 'L''', 'rho*(-l+y-lambda)');
  
  % Pre-solve the KKT system 
  % x = K \ [-f + workLp; b]
  K = [cdr.Q+cdr.rho*cdr.L'*cdr.L cdr.A'; cdr.A zeros(size(cdr.A,1))];
  cdr.mldivide(K, '[-f + workLp; b]', 'workKKT', cdr.KKTSolveMode);  
  cdr.print('x = workKKT(1:size(L,2));\n')

  % Solve prox
  cdr.genProd('workL', cdr.L, 'L', 'x');
  
  cdr.print('q = lambda + workL + l;\n')

  % Evaluate prox functions
  cdr.genProx('q', 'rho', 'y');
  
  % >>>>>>>>>>> COPY >>>>>>>>>>
  
  % Dual update
  lambda = q - y;
  
  % Convergence checks
  if mod(iii, itr_conv) == 0
    s = -rho*L'*(y - prev_y);
    r = lambda - prev_lambda;
    rDual   = s'*s;
    rPrimal = r'*r;
    if rDual < (dualTol*dualTol) && rPrimal < (primalTol*primalTol)
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
rPrimal = sqrt(rPrimal);

%% >>>>>>>>>>> RUN >>>>>>>>>>>
% Must call this to terminate the function (and compute the computation time)
% All addFunction calls must come before this point
cdr.genFunctionEnd();

