classdef normFunc < splitFunc
 %
 % Implements standard norms for splitvars
 %
 % Syntax is the same as Matlab
 %
 % Matrix norms:
 %
 %   norm(X,2) returns the 2-norm of X
 %   norm(X) is the same as norm(X,2)
 %   norm(X,1) returns the 1-norm of X
 %   norm(X,Inf) returns the infinity norm of X
 %   norm(X,'fro') returns the Frobenius norm of X.
 %
 % Vector norms:
 %
 %   norm(V,1) returns the 1-norm of V
 %   norm(V,2) returns the 2-norm of V
 %   norm(V,Inf) returns the largest element of ABS(V).
 %   norm(V,-Inf) returns the smallest element of ABS(V).
 %
 % DCP rule: Input must be affine
 %
 
 properties
  p % type of the norm
  x % affine function whose norm we're taking
 end
 
 methods
  %------------------------------------------------------------------
  % Constructor only defines inputs
  %------------------------------------------------------------------
  function func = normFunc(varargin)
   func = func@splitFunc(varargin{:});
  end
  
  %------------------------------------------------------------------
  % Main initialization function
  %------------------------------------------------------------------
  function func = initialize(func, x, p)
   if nargin < 3, p = 2; end
   
   % Check arguments for validity
   assert(isvector(x), 'splitvar:norm', 'Matrix norms not yet implemented')
   if nargin < 2 || isempty(p), p = 2; end
   
   assert(isa(x, 'splitvar'), 'Can only take the norm of a splitvar')
   if ischar(p)
    assert(strcmp(p,'inf'), 'p must be 1, 2 or inf')
    p = inf;
   else
    assert(isnumeric(p) && isscalar(p), 'p must be 1, 2 or inf')
    assert(p == 1 || p == inf || p == 2, 'p must be 1, 2 or inf')
   end
   
   func.x = x;
   func.p = p;
  end
  
  %------------------------------------------------------------------
  % Return the type of this function
  %------------------------------------------------------------------
  function type = getType(func)
   isMatrix = ~isvector(func.x);
   switch func.p
    case 1
     if isMatrix, type = splitFuncTypes.norm_matrixone;
     else         type = splitFuncTypes.norm_one;
     end
    case 2
     if isMatrix, type = splitFuncTypes.norm_matrixtwo;
     else         type = splitFuncTypes.norm_two;
     end
    case inf
     if isMatrix, type = splitFuncTypes.norm_matrixinf;
     else         type = splitFuncTypes.norm_inf;
     end
    case 'fro'
     type = splitFuncTypes.norm_matrixfro;
   end
  end
  
  %------------------------------------------------------------------
  % Compute the vexity of the function
  %------------------------------------------------------------------
  function vexity = computeVexity(func)
   vexity = 1; % Norms are always convex
  end
  
  
  %------------------------------------------------------------------
  % Convert to optimizable constraint form
  %
  % One normFunc and one splitvar
  %  (Normball or SOC constraint if 2-norm)
  %------------------------------------------------------------------
  function con = flattenAffine(func, w, aff)
   assert(ismember(func.type, splitFuncTypes.vectorNorms, 'legacy'), 'Matrix norms not yet implemented')
   assert(~isempty(aff), 'Constraint norm(x) <= 0 nonsensical')
   con = [];
   
   if aff.isconstant % Linear term is a constant
    %~~~~~~~~~~~~~~~~~~~~~~
    % Norm-ball
    %~~~~~~~~~~~~~~~~~~~~~~
    
    %  w*norm(x) + cnst >= 0
    % convert to
    %  norm(x) <= c = -cnst/w(1)
    
    c = -aff.val / w;
    assert(c > 0, 'Norm-ball constraint has negative size. Will always be infeasible')
    
    con = normBallConSet(func.x, func.p, c);
   else
    if func.p == 2
     %~~~~~~~~~~~~~~~~~~~~~~
     % Second-order cone
     %~~~~~~~~~~~~~~~~~~~~~~
     
     %  w*norm(x) + aff(x) >= 0
     % convert to
     %  norm(x) <= -aff(x)/w
     
     aff = -aff / w;
     
     con = socConSet(func.x, aff);
    end
   end
  end
  
  
  %------------------------------------------------------------------
  % Create epigraph constraint
  %
  %  norm(x,p) <= t
  %
  %------------------------------------------------------------------
  function con = flattenEpi(func)
   assert(ismember(func.type, splitFuncTypes.vectorNorms, 'legacy'), 'Matrix norms not yet implemented')
   
   % Add a cone constraint depending on the type of norm
   switch func.p
    case 1
     % 1-norm epigraph
     % sum(abs(x)) <= t
     con = nonnegConSet(func - sum(abs(func.x)));
     
    case 2
     % Second-order cone
     con = socConSet(func.x, func);
     
    case {'inf', inf}
     % inf-norm epigraph
     %  max(abs(x)) <= t
     %  -t <= x <= t
     con = nonnegConSet([func - func.x;func.x + func]);
     
    otherwise
     error('Unknown norm')
   end
  end
  
  % ----------------------------------------------------------
  % Convert the function to proximal form
  %
  % Returns empty if the proximal form is not closed-form
  %  (if empty, then the generator will use epigraph form instead)
  % ----------------------------------------------------------
  function con = flattenProx(func)
   assert(ismember(func.type, splitFuncTypes.vectorNorms, 'legacy'), 'Matrix norms not yet implemented')

   con = normProxConSet(func.x, func.p);   
  end
  
 end
end
