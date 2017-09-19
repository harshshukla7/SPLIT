classdef weightedSum
 % Represents a function as
 %   sum w(i) * funcs{i}
 
 properties
  funcs % Cell array of functions
  w     % Weights
  types % Function class names (strings)
  x     % Original splitvar
 end
 
 methods
  % Constructor - x must be a scalar splitvar
  function wsum = weightedSum(x)
   assert(isscalar(x), 'getVarTypes only works for scalar functions')
   
   % Get the functions and variables
   ind = find(x.getVars);
   [funcs, functypes] = splitProb.getFunctions(ind);
   indNonlin = ind(functypes ~= splitFuncTypes.affine);
   indLin    = ind(functypes == splitFuncTypes.affine);
   
   funcs = [funcs{:}]';
   
   % Weighting function
   % (assuming that the last func will be the linear part)
   b = x.basis;
   w = b(indNonlin+1);
   
   % Build splitvar to represent linear portion
   bLin = b; bLin(indNonlin+1) = 0;
   if sum(abs(bLin)) > 0
    funcs{end+1} = splitvar(1,1,bLin);
    w = [w;1];
   end
   
   wsum.funcs = funcs;
   wsum.w     = w;
   wsum.x     = x;
  end
  
  %-----------------------------------------------------------------------
  % Simplify the function
  %-----------------------------------------------------------------------
  function w = simplify(w)
   
   %......................................................................
   % 1. Sum of quadratics and linear terms into a single quadratic
   %......................................................................
   w = w.simplifyQuadratic;
  end
  
  %-----------------------------------------------------------------------
  % Flatten the function
  %
  % flat = flatten(w)
  % List of containment constraints
  %-----------------------------------------------------------------------
  function con = flatten(w)
   
   % Default - we don't know anything special, so just add this as an
   % affine constraint + epigraph
   con = nonnegConSet(w.x);
   
   % Search for patterns that we can express as a set containment
   types    = cellfun(@(x)class(x),w.funcs,'uniformoutput',false);
   isAffine = strcmp(types, 'splitvar');
   
   if length(w.funcs) == 2 && any(isAffine) || length(w.funcs) == 1
    % We have one nonlinear function and zero or one affine
    aff = [];
    
    if length(isAffine) == 2 
     if isAffine(1)
      w.funcs = flipud(w.funcs);
      w.w     = flipud(w.w);
     end
     aff      = w.w(2) * w.funcs{2};
    end
    nl       = w.funcs{1};
    nlWeight = w.w(1);
    
    % Ask the nonlinear function if it knows how to parse this form?
    c = nl.flattenAffine(nlWeight, aff);
    if ~isempty(c), con = c; end
   end
  end
 end
 
 % Simplification methods
 methods (Access = private)
  
  %-----------------------------------------------------------------------
  % Combine quadratics and affine terms into a single quadratic
  %-----------------------------------------------------------------------
  function w = simplifyQuadratic(w)
   quadFound = false;
   for i = 1:length(w.w)
    if strcmp(class(w.funcs{i}), 'quadraticFunc')
     quadFound = true; break
    end
   end
   if quadFound == false, return; end
   
   % Combine indQ functions into a single function
   n = splitProb.numVars;
   H = spalloc(n+1,n+1,n);
   indQ = zeros(length(w.w),1)==1;
   for i = 1:length(w.w)
    if strcmp(class(w.funcs{i}), 'quadraticFunc')
     indQ(i) = true;
     H = H + w.w(i)*w.funcs{i}.H;
    elseif strcmp(class(w.funcs{i}), 'splitvar')
     indQ(i) = true;
     B = w.w(i)*w.funcs{i}.basis;
     H(:,1) = H(:,1) + B/2;
     H(1,:) = H(1,:) + B'/2;
    end
   end
   
   if all(indQ)
    % We're left with one quadratic function
    w.funcs  = {quadraticFunc(H)};
    w.w      = 1;
   else
    
    % Remove all linear and quadratic functions
    w.funcs = {w.funcs{~indQ}};
    w.w     = w.w(~indQ);
    
    % Add one new quadratic, which is the sum
    w.funcs{end+1}  = quadraticFunc(H);
    w.w(end+1)      = 1;
   end
  end
 end
end
