classdef minFunc < splitFunc
  % Pointwise minimum of a vector of convex functions
  % 
  % Convex
  
  properties
    x  
  end
  
  
  methods 
    %------------------------------------------------------------------
    % Constructor only defines inputs
    %------------------------------------------------------------------
    function func = minFunc(varargin)
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
      type = splitFuncTypes.min;
    end
    
    %------------------------------------------------------------------
    % Compute the vexity of the function
    %------------------------------------------------------------------
    function vexity = computeVexity(func)
      vexity = 1; % min is always convex
    end
    
    
    %------------------------------------------------------------------
    % Convert to optimizable constraint form
    %------------------------------------------------------------------
    function flatten(func)
      
      if isempty(Y)
        % We're taking the min of X along the rows
        
        % Epigraph variable
        t = splitvar(size(X,2));

        % Compute the min of constant elements directly
        isConX = all(X.isconstant);
        for i = 1:size(X,2)
          if isConX(i)
            t.sub(i) = min(X.sub(:,i));
          else
            t.sub(i) >= X.sub(:,i);
          end
        end
      else
        % We're taking the min of X and Y elementwise

        % min of matrix and scalar
        if isscalar(X) && ~isscalar(Y),     X = X*ones(size(Y));
        elseif ~isscalar(X) && isscalar(Y), Y = Y*ones(size(X));
        end
        
        assert(all(size(X) == size(Y)), 'splitvar:min', 'X and Y must have the same dimension (or one is a scalar)')
        assert(all(Y.vexity >= 0), 'splitvar:min', 'Y must be affine, constant or convex')
        
        % Compute the min of constant elements directly
        isConX = X.isconstant; isConY = Y.isconstant;
        
        % Epigraph variable
        t = splitvar(size(X));
        
        for i = 1:prod(size(X))
          if isConX(i) && isConY(i)
            t = min(X.sub(i).val, Y.sub(i).val)
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
