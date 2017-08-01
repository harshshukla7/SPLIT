%% Define inputs, outputs and parameters

% Parameters available at compile-time and run-time
cdr.addConst('Q','');
cdr.addConst('itr_conv', 'Number of iterations between convergence checks')

% Arguments to the generated function
% cdr.addInput

% Outputs of the generated function
cdr.addOutput('y_',      'Copy of the primal variable');
cdr.addOutput('itr',     'Number of iterations');
cdr.addOutput('tm',      'Computation time');

% Parameters available at compile-time, but not run-time


% Define CP parameters and their defaults
cdr.addParameter('tau',   'tau')
cdr.addParameter('sigma', 'sigma')
cdr.addParameter('gamma', 'gamma')

% define matrices specific for chambolle-pock
Mu        = max(svds([cdr.L; cdr.A]));
cdr.tau   = 1/Mu;
cdr.sigma = 1/(cdr.tau*Mu^2);
cdr.gamma = 0.1*min(eig(cdr.Q));

cdr.addConst('scale_rPrimal', 'Initial primal scaling');
cdr.addConst('scale_rDual',   'Initial dual scaling');

cdr.scale_rPrimal = eye(size(cdr.Q,1));
cdr.scale_rDual   = eye(size(cdr.L,1)+size(cdr.A,1));

% >>>>>>>>>>> IO >>>>>>>>>>>
%
% CPII template
%
% >>>>>>>>>>> HEADER >>>>>>>>>>>


%% >>>>>>>>>>> RUN  >>>>>>>>>>>

% Allocate memory for internal variables
[m,n] = size(cdr.A);
cdr.addVar('x',     n, 1);
cdr.addVar('x_bar', n, 1);
cdr.addVar('nu',    m, 1);
cdr.addVar('p',     size(cdr.L,1), 1);
cdr.print('p_nu  = [p;nu];\n');
cdr.addVar('y',     length(cdr.l),1);

% Compute the parametric solution
cdr.genComputeParametric('par');

% Initial values
cdr.print('tau   = %e;\n', cdr.tau)
cdr.print('sigma = %e;\n', cdr.sigma)
cdr.print('gamma = %e;\n', cdr.gamma)

%% >>>>>>>>>>> COPY >>>>>>>>>>>

rDual   = 0;
rPrimal = 0;
itr     = 0;

for iii = 1:maxItr
  
  % step 0
  p_prev     = p;
  nu_prev    = nu;
  x_prev     = x;
  x_bar_prev = x_bar;
  p_nu_prev  = [p_prev; nu_prev];
  
  % step 1: compute p and nu
  p  = p  + sigma*(L*x_bar + l);
  nu = nu + sigma*(A*x_bar - b);
  
  %% >>>>>>>>>>> RUN >>>>>>>>>>>
  % step 2: compute projection
  cdr.genProx('p', 'sigma', 'y');
  %% >>>>>>>>>>> COPY >>>>>>>>>>>
  
  % step 3: compute x
  % this formula represents the conventional procedure
  % x = (dat.Q + eye(size(dat.Q, 1))*(1/tau)) \ ((1/tau)*x - dat.f' - (K'*p_nu));
  
  % exploit structure of matrix Q - simpler C implementation
  % pull out elements from the diagonal and update x  
  diagQ    = diag(Q);
  diagItau = diag(eye(size(Q, 1))*(1/tau));
  inverse  = 1./(diagQ + diagItau);
  x = inverse .* ((1/tau)*x - f - ([L' A']*[p;nu]));
  
  % step 4: update stepsize
  theta = 1 / (sqrt(1+2*gamma*tau));
  tau   = theta*tau;
  sigma = sigma / theta;
  
  % step 5: update x_bar
  x_bar = x + theta*(x - x_prev);
  
  % step 6: convergence check
  if mod(iii, itr_conv) == 0
    invSP = inv(scale_rPrimal);
    invSD = inv(scale_rDual);
    K = [L; A];
    
    r  = (1/tau)*(invSP*(x_prev-x)) - (K*scale_rPrimal)'*(p_nu_prev - p_nu);
    s  = (1/sigma)*scale_rDual'*(p_nu_prev-p_nu) - invSD*K*(x - x_bar);
    
    rPrimal = r'*r;
    rDual   = s'*s;
    
    if rDual < dualTol^2 && rPrimal < primalTol^2
      break
    end
  end
end

% Set the output
itr     = iii;
x_      = x;
y_      = y;

%% >>>>>>>>>>> RUN >>>>>>>>>>>
cdr.genFunctionEnd();
