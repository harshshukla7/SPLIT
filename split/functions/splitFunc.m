classdef splitFunc < splitvar
 %
 % splitFunc
 %
 % Represents a nonlinear function of splitvars via an epigraph
 %
 % Abstract class. Specific functions are derived from this
 %
 
 properties
  type   % splitFuncTypes
 end
 
 methods
  function func = splitFunc(varargin)
   if nargin == 0
    args = {1,0,0};
   else
    args = {};
   end
   
   func = func@splitvar(args{:});
   
   if nargin > 0
    % Let the specific function type initialize itself
    func = func.initialize(varargin{:});
    
    % Get the function type
    func.type = func.getType;
    
    % Compute vexity
    vexity = func.computeVexity;
    
    % Validate disciplined convex programming ruleset
    assert(func.validateDCP, 'Function does not satisfy the DCP rules')
    
    % Register this function with splitProb
    splitProb.registerFunction(func, vexity);
   end
  end
  
  % Validate that the function type(x,y) satisfies the DCP rules
  function isValid = validateDCP(func)
   %%% TODO
   isValid = true;
  end
 end
 
 methods
  % ----------------------------------------------------------
  % Convert the function + affine to
  %   linear constraints + set containment
  %
  % Overload in sub-functions if this is possible
  % ----------------------------------------------------------
  function con = flattenAffine(func, w, aff)
   con = [];
  end
  
  % ----------------------------------------------------------
  % Convert the function to epgraph form
  %
  % Must be overloaded in sub-functions
  % ----------------------------------------------------------
  function con = flattenEpi(func)
   error('flattenEpi has not yet been implemented for functions of type %s', char(func.type))
  end
  
  % ----------------------------------------------------------
  % Convert the function to proximal form
  %
  % Returns empty if the proximal form is not closed-form
  %  (if empty, then the generator will use epigraph form instead)
  % ----------------------------------------------------------
  function con = flattenProx(func)
   con = [];
  end
 end
 
 methods (Abstract)
  initialize(func)
  computeVexity(func)
  getType(func)
 end
end
