classdef absFunc < splitFunc
 %
 % abs(x)
 %
 
 properties
  x % affine function whose norm we're taking
 end
 
 methods
  %------------------------------------------------------------------
  % Constructor only defines inputs
  %------------------------------------------------------------------
  function func = absFunc(varargin)
   func = func@splitFunc(varargin{:});
  end
  
  %------------------------------------------------------------------
  % Main initialization function
  %------------------------------------------------------------------
  function func = initialize(func, x)
   assert(isa(x, 'splitvar'), 'Can only take the absolute value of a splitvar')
   assert(isscalar(x), 'Can only take the absolute value of a scalar')
   
   func.x    = x;
  end
  
  %------------------------------------------------------------------
  % Return the type of this function
  %------------------------------------------------------------------
  function type = getType(func)
   type = splitFuncTypes.abs;
  end
  
  %------------------------------------------------------------------
  % Compute the vexity of the function
  %------------------------------------------------------------------
  function vexity = computeVexity(func)
   vexity = 1; % Abs is always convex
  end
  
  % ----------------------------------------------------------
  % Convert the function + affine to
  %   linear constraints + set containment
  %
  %   w|x| + aff >= 0
  %
  % Convert to two linear constraints
  %
  %   w*x + aff >= 0
  %  -w*x + aff >= 0
  %
  % ----------------------------------------------------------
  function con = flattenAffine(func, w, aff)
   con = nonnegConSet([w*func + aff;-w*func + aff]);
  end

  %------------------------------------------------------------------
  % Create epigraph constraint
  %
  %  |x| <= t 
  %
  % x is a scalar
  %
  %------------------------------------------------------------------
  function con = flattenEpi(func)
   % -t <= x <= t
   con = nonnegConSet([func - func.x;func.x + func]);
  end
  
  
  %       t = splitvar(size(x));
  %       isCon = x.isconstant;
  %
  %       if any(isCon(:))
  %         % Directly compute the abs of constants
  %         t = asgn(t, abs(x.sub(isCon).val), isCon);
  %       end
  %
  %       % Option 1 : Implement using lorentz cone
  %       %       for i = 1:prod(size(t))
  %       %         if isCon(i), continue; end
  %       %         inSecondOrderCone(x.sub(i), t.sub(i));
  %       %       end
  %
  %       % Option 2 : Implement using upper / lower bounds
  %       -t.sub(~isCon) <= x.sub(~isCon) <= t.sub(~isCon);
  %
  %       splitProb.setVexity(t.sub(~isCon).getVarInd, 1);
  %     end
 end
end
