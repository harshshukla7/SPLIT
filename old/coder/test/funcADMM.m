function [x_, y_, lambda_, itr, tm, rDual] = funcADMM(par, f, b, l, A, L, Q, pB, pF, pL, dualTol, primalTol, maxItr, rho, itr_conv, iK_514, iK_402)%#codegen
%
% ADMM template
%
% Add any custom comments here
%
% This whole block must be comments (no spaces), so that the automatically
% generated comments appear when 'help xxx' is typed at the command line
%
%
% Template file:             /Users/cnjones/svn/split/split/coder/@splitCoder/templates/template_admm.m
% Generated matlab function: funcADMM
% Generated mex function:    funcADMM_mex
%
% Problem:
%  3 variables
%  2 parameters
%  1 proximal functions:
%    - nonNegative
%
% Function prototype:
%  [x_, y_, lambda_, itr, tm, rDual] = funcADMM(par)
%  [x_, y_, lambda_, itr, tm, rDual] = funcADMM_mex(par)
%
% Arguments to the generated function
%                 par %  2 x 1    Parameter variable
%
% Outputs of the generated function
%                  x_ %  ? x ?    Optimizer
%                  y_ %  ? x ?    Copy of the primal variable
%             lambda_ %  ? x ?    Dual variable
%                 itr %  ? x ?    Number of iterations
%                  tm %  ? x ?    Computation time
%               rDual %  ? x ?    Dual residual
%
% Parameters available at compile-time and run-time
%                   f %  3 x 1    
%                   b %  1 x 1    
%                   l % 20 x 1    
%                   A %  1 x 3    
%                   L % 20 x 3    
%                   Q %  3 x 3    
%                  pB %  1 x 2    
%                  pF %  3 x 2    
%                  pL % 20 x 2    
%             dualTol %  1 x 1    Dual tolerance
%           primalTol %  1 x 1    Primal tolerance
%              maxItr %  1 x 1    Maximum number of iterations
%                 rho %  1 x 1    Penalty parameter
%            itr_conv %  1 x 1    Number of iterations between convergence checks
%              iK_514 %  4 x 4    Pre-computed matrix inverse
%              iK_402 %  4 x 4    Pre-computed matrix inverse
%
% Parameters available at compile-time, but not run-time
%        KKTSolveMode %  1 x 10   'LU','Invert KKT'
%  MatVec_SparseLimit %  1 x 1    Percent sparsity below which matrix-vector products use code-gen
%

if coder.target('MEX'), coder.ceval('split_tic'); else tic; end

persistent x y lambda workLp workKKT workL 
if isempty(x)
	x = zeros(size(L,2), 1);
	y = zeros(20, 1);
	lambda = zeros(20, 1);
	workLp = zeros(3, 1);
	workKKT = zeros(4, 1);
	workL = zeros(20, 1);
end


f = pF*par + f;
rDual   = 0;
rPrimal = 0;
itr     = 0;

for iii = 1:maxItr
  prev_lambda = lambda;
  prev_y      = y;
  
  % Solve x+
workLp = L' * rho*(-l+y-lambda);
workKKT = iK_402*[-f + workLp; b];
x = workKKT(1:size(L,2));
workL = L * x;
q = lambda + workL + l;
y = proxFunc(q, rho, y);
  
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

tm = 0.0;
if coder.target('MEX')
  tm = coder.ceval('split_toc');
else
  tm = toc;
end
end


function y = proxFunc(q, rho, y)
y(1:20) = proj_positive(q(1:20));
end
