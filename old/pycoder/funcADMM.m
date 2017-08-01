function [x_, y_, lambda_, itr, tm, rDual] = funcADMM(par, f, b, l, A, L, Q, pB, pF, pL, dualTol, primalTol, maxItr, rho, itr_conv, iK_240)%#codegen
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
%  19 variables
%  2 parameters
%  2 proximal functions:
%    - nonNegative
%    - secondOrderCone
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
%                   f % 19 x 1    
%                   b % 10 x 1    
%                   l % 41 x 1    
%                   A % 10 x 19   
%                   L % 41 x 19   
%                   Q % 19 x 19   
%                  pB % 10 x 2    
%                  pF % 19 x 2    
%                  pL % 41 x 2    
%             dualTol %  1 x 1    Dual tolerance
%           primalTol %  1 x 1    Primal tolerance
%              maxItr %  1 x 1    Maximum number of iterations
%                 rho %  1 x 1    Penalty parameter
%            itr_conv %  1 x 1    Number of iterations between convergence checks
%              iK_240 % 29 x 29   Pre-computed matrix inverse
%
% Parameters available at compile-time, but not run-time
%        KKTSolveMode %  1 x 10   'LU','Invert KKT'
%  MatVec_SparseLimit %  1 x 1    Percent sparsity below which matrix-vector products use code-gen
%

if coder.target('MEX'), coder.ceval('split_tic'); else tic; end

persistent x y lambda workLp workKKT workL 
if isempty(x)
	x = zeros(size(L,2), 1);
	y = zeros(41, 1);
	lambda = zeros(41, 1);
	workLp = zeros(19, 1);
	workKKT = zeros(29, 1);
	workL = zeros(41, 1);
end


b = pB*par + b;
rDual   = 0;
rPrimal = 0;
itr     = 0;

for iii = 1:maxItr
  prev_lambda = lambda;
  prev_y      = y;
  
  % Solve x+
workLp = workLp_func_184(rho*(-l+y-lambda), workLp);
workKKT = iK_240*[-f + workLp; b];
x = workKKT(1:size(L,2));
workL = workL_func_418(x, workL);
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


function y = workLp_func_184(x, y)
y(1) =  + x(1) - x(11);
y(2) =  + x(2) - x(12);
y(3) =  + x(3) - x(13);
y(4) =  + x(4) - x(14);
y(5) =  + x(5) - x(15) + x(40);
y(6) =  + x(6) - x(16) + x(41);
y(7) =  + x(7) - x(17) - x(30) + x(32);
y(8) =  + x(8) - x(18) - x(31) + x(33);
y(9) =  + x(9) - x(19) - x(35) + x(36);
y(10) =  + x(10) - x(20) - x(37) + x(38);
y(11) =  + x(21) - x(25);
y(12) =  + x(22) - x(26);
y(13) =  + x(23) - x(27);
y(14) =  + x(24) - x(28);
y(15) =  - x(29) + x(39);
y(16) =  - x(29) + x(30) + x(31) + x(32) + x(33);
y(17) =  - x(29) + x(34);
y(18) =  - x(34) + x(35) + x(36);
y(19) =  - x(34) + x(37) + x(38);
end


function y = workL_func_418(x, y)
y(1) = +x(1);
y(2) = +x(2);
y(3) = +x(3);
y(4) = +x(4);
y(5) = +x(5);
y(6) = +x(6);
y(7) = +x(7);
y(8) = +x(8);
y(9) = +x(9);
y(10) = +x(10);
y(11) = -x(1);
y(12) = -x(2);
y(13) = -x(3);
y(14) = -x(4);
y(15) = -x(5);
y(16) = -x(6);
y(17) = -x(7);
y(18) = -x(8);
y(19) = -x(9);
y(20) = -x(10);
y(21) = +x(11);
y(22) = +x(12);
y(23) = +x(13);
y(24) = +x(14);
y(25) = -x(11);
y(26) = -x(12);
y(27) = -x(13);
y(28) = -x(14);
y(29) =  - x(15) - x(16) - x(17);
y(30) =  - x(7) + x(16);
y(31) =  - x(8) + x(16);
y(32) =  + x(7) + x(16);
y(33) =  + x(8) + x(16);
y(34) =  + x(17) - x(18) - x(19);
y(35) =  - x(9) + x(18);
y(36) =  + x(9) + x(18);
y(37) =  - x(10) + x(19);
y(38) =  + x(10) + x(19);
y(39) = +x(15);
y(40) = +x(5);
y(41) = +x(6);
end


function y = proxFunc(q, rho, y)
y(1:38) = proj_positive(q(1:38));
y(39:41) = proj_secondOrderCone(q(39:41));
end
