classdef maxFunc < splitFunc
  % Pointwise maximum of a vector of convex functions
  % 
  % Convex
  
  properties
    x  
  end
  
  
  methods 
    %------------------------------------------------------------------
    % Constructor only defines inputs
    %------------------------------------------------------------------
    function func = maxFunc(varargin)
      func = func@splitFunc(varargin{:});
    end
    
    %------------------------------------------------------------------
    % Main initialization function
    %------------------------------------------------------------------
    function func = initialize(func, x)
      assert(isa(x, 'splitvar'), 'x must be a splitvar')
      assert(all(x.vexity >= 0), 'X must be affine, constant or convex')
      
      func.x = x;
    end

    %------------------------------------------------------------------
    % Return the type of this function
    %------------------------------------------------------------------    
    function type = getType(func)
      type = splitFuncTypes.max;
    end
    
    %------------------------------------------------------------------
    % Compute the vexity of the function
    %------------------------------------------------------------------
    function vexity = computeVexity(func)
      vexity = 1; % max is always convex
    end
    
    
    %------------------------------------------------------------------
    % Convert to optimizable constraint form
    %------------------------------------------------------------------
    function flatten(func)
      
      if isempty(Y)
        % We're taking the max of X along the rows
        
        % Epigraph variable
        t = splitvar(size(X,2));

        % Compute the max of constant elements directly
        isConX = all(X.isconstant);
        for i = 1:size(X,2)
          if isConX(i)
            t.sub(i) = max(X.sub(:,i));
          else
            t.sub(i) >= X.sub(:,i);
          end
        end
      else
        % We're taking the max of X and Y elementwise

        % Max of matrix and scalar
        if isscalar(X) && ~isscalar(Y),     X = X*ones(size(Y));
        elseif ~isscalar(X) && isscalar(Y), Y = Y*ones(size(X));
        end
        
        assert(all(size(X) == size(Y)), 'splitvar:max', 'X and Y must have the same dimension (or one is a scalar)')
        assert(all(Y.vexity >= 0), 'splitvar:max', 'Y must be affine, constant or convex')
        
        % Compute the max of constant elements directly
        isConX = X.isconstant; isConY = Y.isconstant;
        
        % Epigraph variable
        t = splitvar(size(X));
        
        for i = 1:prod(size(X))
          if isConX(i) && isConY(i)
            t = max(X.sub(i).val, Y.sub(i).val)
          else
            t >= X.sub(i);
            t >= Y.sub(i);
          end
        end
      end
      
      % Define this as an epigraph variable by setting its convex flag
      splitProb.setVexity(t.getVarInd, 1);
    end
  end
end
