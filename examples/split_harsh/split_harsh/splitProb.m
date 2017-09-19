classdef splitProb < handle
 %
 % Contains all data for the problem being created
 %
 
 properties
  nextVarID    = 1;
  
  % List of structs representing the variables containing:
  %  - value      Current value (double)
  %  - vexity     -1, 0, 1 concave / affine / convex
  %  - parameter  true if this is a parameter
  %  - func       splitvar / splitFunc
  %  - type       splitFuncTypes
  %  - epi        0 = not used, -1 = to be processed, 1 = processed
  vars         = [];
  inequalities = [];
  equalities   = [];
  objective    = [];
 end
 
 % Public methods
 methods
  % Add a constraint to the constraint set
  %
  %  a op b
  % where a and b are splitvars and op is one of
  %  'le', 'lt', 'ge', 'gt', 'eq'
  %
  % -- or --
  %
  %  a 'in' b
  % where b is a splitSet (not implemented yet)
  %
  % -- or --
  %
  %  A op B
  % where A and B are symmetric matrices and
  % op is 'succeq', 'preceq'
  %
  function addConstraint(sp, a, b, op)
   if isnumeric(a), a = splitvar(size(a), [], a(:)'); end
   if isnumeric(b), b = splitvar(size(b), [], b(:)'); end
   
   if isscalar(a) && ~isscalar(b), a = ones(size(b))*a; end
   if isscalar(b) && ~isscalar(a), b = ones(size(a))*b; end
   
   assert(all(a.size == b.size), 'Functions must be same size')
   
   vexA = a.vexity;
   vexB = b.vexity;
   
   switch op
    case {'lt', 'le'}
     op = 'lt';
     
     % a must be convex and b concave
     assert(all(vexA(:) >= 0), 'LHS must be convex')
     assert(all(vexB(:) <= 0), 'RHS must be concave')
     
     func = b-a;
     if isempty(sp.inequalities), sp.inequalities = func(:);
     else                         sp.inequalities = [sp.inequalities;func(:)];
     end
     
    case {'gt', 'ge'}
     op = 'gt';
     
     % a must be convex and b concave
     assert(all(vexA(:) <= 0), 'LHS must be concave')
     assert(all(vexB(:) >= 0), 'RHS must be convex')
     
     func = a-b;
     if isempty(sp.inequalities), sp.inequalities = func(:);
     else                         sp.inequalities = [sp.inequalities;func(:)];
     end
     
    case 'eq'
     
     % a and b must be affine
     assert(all(vexA(:) == 0), 'LHS must be affine')
     assert(all(vexB(:) == 0), 'RHS must be affine')
     
     func = a-b;
     if isempty(sp.equalities), sp.equalities = func(:);
     else                       sp.equalities = [sp.equalities;func(:)];
     end
     
    case 'succeq'
     
     error('PSD constraints not implemented yet')
     
    case 'preceq'
     
     error('PSD constraints not implemented yet')
     
    case 'in'
     error('Operator ''in'' not implemented yet')
     
    otherwise
     error('Unknown operation %s', op)
   end
  end
 end
 
 % Make the constructor private, so that we get a singleton
 methods (Access=private)
  function this = splitProb
  end
 end
 
 methods(Static)
  % Get the object via the static method instance
  function obj = instance(~)
   persistent uniqueInstance
   if isempty(uniqueInstance) || nargin > 0
    obj = splitProb();
    uniqueInstance = obj;
    
    % Print tasks to do
    todo
   else
    obj = uniqueInstance;
   end
  end
  
  % Clear the problem
  function clearProblem
   splitProb.instance(true);
  end
  
  % Return total number of variables
  function nvars = numVars
   %    nvars = splitProb.instance.nextVarID - 1;
   nvars = length(splitProb.instance.vars);
  end
  
  % Get unique variable indices
  function varID = getNewVarID(dim)
   if nargin < 1, dim = 1; end
   dim = dim(:);
   assert(length(dim) <= 2, 'Variable can only have up to two dimensions')
   if length(dim) == 1, dim = [dim 1]; end
   n = prod(dim);
   
   sp = splitProb.instance;
   
   varID = reshape(sp.nextVarID+[0:n-1], dim(1), dim(2));
   sp.nextVarID  = sp.nextVarID + n;
   
   % Generate new affine variables
   [dat(1:n).func]      = deal([]);
   [dat(1:n).value]     = deal(0);
   [dat(1:n).vexity]    = deal(0);
   [dat(1:n).parameter] = deal(0);
   [dat(1:n).type]      = deal(splitFuncTypes.affine);
   [dat(1:n).epi]       = deal(0);
   
   sp.vars = [sp.vars;dat(:)];
  end
  
  %--------------------------------------------------------------------
  % Flag a variable as an epigraph variable
  % Input: splitFunc
  %--------------------------------------------------------------------
  function registerFunction(func, vexity)
   sp = splitProb.instance;
   ind = func.varNum_;
   sp.vars(ind).vexity = vexity;
   sp.vars(ind).func   = {func};
   sp.vars(ind).type   = func.type;
   sp.vars(ind).epi    = 0;
  end
  
  function vex = getVexity(ind)
   %--------------------------------------------------------------------
   % get.vexity
   % Get the vexity of a subset of the variables
   %
   %--------------------------------------------------------------------
   sp  = splitProb.instance;
   vex = [sp.vars.vexity]';
   if nargin > 0, vex = vex(ind); end
  end
  
  function add(a, b, op)
   %--------------------------------------------------------------------
   % Add a constraint to the stack
   %--------------------------------------------------------------------
   splitProb.instance.addConstraint(a,b,op);
  end
  
  function [funcs,types,isNonlinear] = getFunctions(ind)
   %--------------------------------------------------------------------
   % Get the functions for the given indices
   % Returns a cell-array
   %--------------------------------------------------------------------
   sp = splitProb.instance;
   if nargin == 0, ind = [1:length(sp.vars)]; end
   sp = splitProb.instance;
   types       = [sp.vars(ind).type]';
   funcs       = {sp.vars(ind).func}';
   isNonlinear = ind(find(types ~= splitFuncTypes.affine));
  end
  
  function vals = varValues
   %--------------------------------------------------------------------
   % Set field of vars structure to given val vector
   %--------------------------------------------------------------------
   sp = splitProb.instance;
   vals = [sp.vars.value]';
  end
  
  function setObjective(x)
   %--------------------------------------------------------------------
   % Set the function to be minimizes
   %--------------------------------------------------------------------
   assert(isa(x, 'splitvar') && isscalar(x), 'splitProb:argchk', 'Function to be minimized must be a scalar')
   sp = splitProb.instance;
   sp.objective = x;
  end
  
  function setVarValue(ind, values)
   %--------------------------------------------------------------------
   % Set the current values of the variables indexed by ind
   %--------------------------------------------------------------------
   sp = splitProb.instance;
   assert(numel(ind) == numel(values) || isscalar(values), 'Number of values specified must match number of variables to set');
   assert(max(ind) <= splitProb.numVars, 'Attempt to set the value of a non-existent variable')
   sp.setVars('value', ind, values);
  end
  
  function vals = getVarValue(ind)
   %--------------------------------------------------------------------
   % Get the current values of the variables indexed by ind
   %--------------------------------------------------------------------
   sp   = splitProb.instance;
   vals = [sp.vars(ind).value]';
  end
  
  function setSolution(prob, sol)
   %--------------------------------------------------------------------
   % Set the current value of the variables to the provided solution
   %--------------------------------------------------------------------
   sp   = splitProb.instance;
   % Zero out all variables
   sp.setVars('value', [1:sp.numVars],zeros(sp.numVars,1));
   % Set the solution
   sp.setVars('value', prob.vars, sol.x);
  end
  
  function setVarToParameter(ind)
   %--------------------------------------------------------------------
   % Flag variables as parameters
   %--------------------------------------------------------------------
   sp = splitProb.instance;
   for i = 1:length(ind)
    sp.vars(ind(i)).parameter = 1;
   end
  end

 end
 
 methods
  function setVars(sp, field, ind, val)
   %--------------------------------------------------------------------
   % Set field of vars structure to given val vector
   %--------------------------------------------------------------------
   v = num2cell(val);
   [sp.vars(ind).(field)] = v{:};
  end
 end
 
 
 %% Main problem generation methods
 methods (Static)
  %%
  function prob = genProblem
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Generate the problem data from a flattened description
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %
   %  Output format (assuming no parameters):
   %
   %   min  0.5*x'*Q*x + f*x + c + sum l_i (L*x + l)
   %   s.t.  A*x = b
   %
   %
   %  Output format (with parameters y):
   %
   %   min  0.5*x'*Q*x + (y'*pF + f)*x + (y'*pQ*y + pf'*y + c)
   %                          + sum l_i (L*x + pL*y + l)
   %   s.t.  A*x = b + pB*y
   %
   %
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   sp = splitProb.instance;
   dat.objective    = sp.objective;
   dat.inequalities = sp.inequalities;
   dat.equalities   = sp.equalities;
   
   % Convert from hierarchical to weighted-sum form
   flat = sp.flatten(dat);   
   % Generate problem data
   prob = sp.genProbData(flat);
   
   prob.input.dat  = dat;
   prob.input.flat = flat;
  end

 end
 
 methods (Static)
  function flat = flatten(dat)
   %
   % Convert the hierarchical structure used for construction to a weighted
   % sum of functions
   %
   
   sp = splitProb.instance;
   
   %.....................................................................
   % Process the objective function
   %.....................................................................
   assert(~isempty(dat.objective), 'Objective not specified')
   obj = weightedSum(dat.objective).simplify;
   % Seperate in a quadratic/linear term + nonlinear terms
   ind = cellfun(@(x)x.type == splitFuncTypes.quadratic || x.type == splitFuncTypes.affine, obj.funcs);
   assert(sum(ind) <= 1, 'There should not be multiple linear or quadratic terms in objective after simplification')
   % Nonlinear and non-quadratic terms in the objective function
   nlObj   = {obj.funcs{~ind}}; w = [obj.w(~ind)];
   % Quadratic / affine portion of the objective function
   quadObj = obj.funcs{ind};
   
   % Objective is now: quad + sum w(i)*nlObj(i)
   % For each nlObj, convert to a prox form if possible, otherwise convert
   % to an epigraph form

   con = {};
   for i = 1:length(nlObj)
    tmp = nlObj{i}.flattenProx;
    if ~isempty(tmp)
     tmp.weight = w(i);
     w(i) = 0; % This one has now been delt with, so we're not going to include it in the epigraph list     
     con{end+1} = tmp;
    else
     % It doesn't have a nice prox form - so turn it into an epigraph
     sp.vars(nlObj.getVars(true)) = -1;
    end
   end
   
   % Create linear basis for the nonlinear terms in epigraph form
   alpha = spalloc(sp.numVars+1,1,length(obj.funcs));
   ind   = [cellfun(@(x)x.varNum_, nlObj)];
   alpha(ind+1) = w;
   
   % Reduce to a quadraticFunc
   flat.obj = quadObj.addLinear(alpha);
   
   % Flag all the non-linear terms in the object as epigraph variables
%    for i = 1:length(ind), sp.vars(ind(i)).epi = -1; end
   
   %.....................................................................
   % Grab all the affine terms
   %.....................................................................
   aff = [];
   if ~isempty(dat.inequalities)
    affInd  = dat.inequalities.isaffine;
    if any(affInd), aff     = dat.inequalities(affInd); end

    if any(~affInd)
     % Break the remaining non-linear constraints into weighted sums and
     % then search for known set-containment patterns
     nonlin  = dat.inequalities(~affInd);
     for i = 1:length(nonlin)
      con{end+1} = weightedSum(nonlin(i)).simplify.flatten;
     end
     
     % Add any nonlinear functions to the epigraph list
     e = unique(cell2mat(cellfun(@(x)x.epi',con,'uniformoutput',false)));
     [~,~,epiInd] = splitProb.getFunctions(e);
     [sp.vars(epiInd).epi] = deal(-1); % Flag nonlinear variables for epigraph conversion
     
     % Remove empty constraint sets (they were only there to grab the epigraph constraints)
     con = {con{cellfun(@(x)x.type ~= conSetTypes.noConstraint, con)}};
     % Grab any simple linear inequalities
     ind = find(cellfun(@(x)x.type == conSetTypes.nonNegative, con));
     for i = 1:length(ind)
      aff = [aff;con{ind(i)}.aff(:)];
     end
     con = {con{cellfun(@(x)x.type ~= conSetTypes.nonNegative, con)}};
    end
   end
   
   %.....................................................................
   % Process all used epigraph variables into set-containments
   %.....................................................................
   while 1
    % Grab an epigraph variable marked to be processed
    ind = find([sp.vars.epi]==-1,1);
    if isempty(ind), break; end
    
    func = sp.vars(ind).func{1};
    c  = func.flattenEpi;
    assert(c.type ~= conSetTypes.noConstraint, 'Function %s was not able to flatten its epigraph', class(func))
    sp.vars(ind).epi = 1; % Mark this variable as processed
    
    % Grab linear constraints
    if c.type == conSetTypes.nonNegative
     aff = [aff;c.aff];
    else
     con{end+1} = c;
    end
       
    % Mark the epi variables for processing
    [~,~,epiInd] = splitProb.getFunctions(c.epi);
    for i = 1:length(epiInd)
     if sp.vars(epiInd(i)).epi == 0
      sp.vars(epiInd(i)).epi = -1;
     end
    end
   end
   
   %.....................................................................
   % Finalize the flattened problem
   %.....................................................................
   flat.con = con;
   flat.aff = aff;
   flat.eq  = sp.equalities;
  end

  function prob = genProbData(flat)
   %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   % Compute the optimization problem data from a flattened description
   %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
   %.....................................................................
   % Generate equality constraints
   %.....................................................................
   if ~isempty(flat.eq)
    B = flat.eq.basis;
    dat.A = B(2:end,:)';
    dat.b = -B(1,:)';
   else
    dat.A = zeros(0,splitProb.numVars);
    dat.b = zeros(0,1);
   end
   
   %.....................................................................
   % Generate cost function
   %.....................................................................
   if flat.obj.type == splitFuncTypes.quadratic
    H = flat.obj.H;
    dat.Q = 2*H(2:end,2:end);
    dat.f = 2*H(2:end,1)';
    dat.c = H(1,1);
   else
    B = flat.obj.basis;
    splitProb.numVars;
    dat.Q = spalloc(n,n,0);
    dat.f = B(2:end)';
    dat.c = B(1);
   end
   
   %.....................................................................
   % Generate inequality constraints
   %.....................................................................
   prox = [];
   if length(flat.aff) > 0
    prox = [prox;nonnegConSet(flat.aff).proxData];
   end
   
   %.....................................................................
   % Generate proximal functions
   %.....................................................................
   for i = 1:length(flat.con)
    prox = [prox;flat.con{i}.proxData];
   end
   
   %.....................................................................
   % Remove unused variables
   %.....................................................................
   used = zeros(1,splitProb.numVars)==1;
   % Equality constraints and cost
   used = any([used;any(dat.A,1);any(dat.f,1);any(dat.Q)]);
   % Constraints
   for i = 1:length(prox)
    used = any([used;any(prox(i).L)]);
   end
   unused = ~used;
   
   % Removed unused variables
   dat.A(:,unused) = [];
   dat.f(:,unused) = [];
   dat.Q = dat.Q(used,used);
   for i = 1:length(prox)
    prox(i).L(:,unused) = [];
   end   
   
   prob.vars = find(used); % Mark used variables for recovery
   prob.dat  = dat;
   prob.prox = prox;   

   % Convert to parametric problem (if any parametric variables)
   prob = splitProb.processParametericProb(prob);
  end
  
  function prob = processParametericProb(prob)
   %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   % Pull out all parametric variables
   %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   sp = splitProb.instance;
   
   % Get list of parametric variables
   prob.param = find([sp.vars.parameter] == 1);

   if isempty(prob.param), return; end % No parametric variables

   % Find each parametric variable in the existing list in the right order,
   % so that the columns of the existing matrices correspond to the right
   % parameter
   indParam = zeros(length(prob.param),1);
   for i = 1:length(prob.param)
    indParam(i) = find(prob.vars == prob.param(i));
   end
   varsTmp = prob.vars;
   prob.vars = setdiff(prob.vars, prob.param);   
   indVars  = zeros(length(prob.vars)-length(prob.param),1);
   for i = 1:length(prob.vars)
    indVars(i) = find(varsTmp == prob.vars(i));
   end
      
   % Process the linear equality constraints
   prob.dat.pB = -prob.dat.A(:,indParam);
   prob.dat.A(:,indParam) = [];
   
   % Process the quadratic objective function
   %  0.5*x'*Q*x + (y'*pF + f)*x + (0.5*y'*pQ*y + pf'*y + c)
   %                               ^^^^^^^^^^^^^^^^^^^^^^^^
   %                                 Constant
   prob.dat.pQ =   prob.dat.Q(indParam,indParam);
   prob.dat.pf =   prob.dat.f(indParam);
   
   prob.dat.pF =   prob.dat.Q(indVars,indParam)';
   prob.dat.f  =   prob.dat.f(indVars);
   prob.dat.Q  =   prob.dat.Q(indVars,indVars);
   
   % Process the prox operators
   for i = 1:length(prob.prox)
    prob.prox(i).pL = prob.prox(i).L(:,indParam);
    prob.prox(i).L(:,indParam) = [];
   end
  end
 end
end
