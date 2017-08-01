classdef splitvar
% TODO:
%  - add a disp function. Include checks for mixed parameter / splitvars
  
  
  properties
    size_   = [0 0];
    basis_  = [];
    varNum_ = []; % Global index in splitProb for this variable
  end
  
  methods
    % Creates a new variable vector of given size
    function obj = splitvar(n, m, basis)
      if nargin < 3, basis = NaN; end
      if nargin < 2, m = 1; end
      if nargin < 1, n = 1; end

      assert(~isempty(basis), 'Empty basis specified - this will incorrectly create a new variable')
      
      if numel(n) == 2, sz = n; 
      else              sz = [n m];
      end

      assert(sum(sz) > 0, 'Cannot create a splitvar of size 0x0')
      
      obj.size_  = sz;
      if ~isnan(basis)
        if isnumeric(basis)
          assert(size(basis,2) == prod(sz), 'splitvar:sizechk', 'Matrix dimensions must agree')
          obj.basis_ = basis;
        elseif ischar(basis) && strcmp(basis, 'symmetric')
          % Create a symmetric matrix
          assert(sz(1) == sz(2), 'Symmetric matrices must be square')
          n = sz(1);
          tmp = splitvar(n*(n+1)/2);
          ind = symMat([1:n*(n+1)/2]); ind = ind(:);
          obj.size_ = sz;
          obj.basis_ = tmp.basis_(:,ind);
        else
          error('Unknown type for the third argument of splitvar constructor')
        end
      else
        % Create new variables
        varNum      = splitProb.getNewVarID(obj.size_);
        obj.varNum_ = varNum;
        obj.basis_  = spalloc(splitProb.numVars+1,prod(obj.size_),ceil(0.1*prod(obj.size_)));
        varNum      = varNum(:);
        obj.basis_  = sparse(varNum+1,1:length(varNum),1);
      end
    end
        
    % Overload basic matrix requests
    function n  = numel(varargin), n = 1; end
    function n  = numelx(obj)
      if isempty(obj), n = 0; else n  = prod(obj.size_); end
    end
    function [sz,n] = size(obj, dim)
      if nargin<2, dim = [1 2]; end
      if isempty(obj)
        sz = builtin('size',obj);
      else
        sz = obj.size_(dim);
      end
      
      if nargout > 1
        n = sz(2); sz = sz(1);
      end
    end

    function l = length(obj)
      if isempty(obj)
        l = 0;
      else
        l = max(obj.size);
      end
    end
    
    % Get the basis, but lift it to the current size of the defined variables first
    function b = basis(obj)
      b = obj.basis_;
      if size(b,1) < splitProb.numVars+1
        b = [b;spalloc(splitProb.numVars+1-size(b,1),size(b,2),3*size(b,2))];
      end
    end
  end
  
  %--------------------------------------------------
  %--------------------------------------------------
  % Basic arithmetic operators
  %--------------------------------------------------
  %--------------------------------------------------
  methods
    
    %--------------------------------------------------
    % Selecting sub-matrices
    %
    % Shorthand for subsref for use within member functions
    %--------------------------------------------------
    function mat = sub(x,I,J)
      S.type = '()';
      if nargin < 3, S.subs = {I};
      else           S.subs = {I,J};
      end
      mat = subsref(x,S);
    end
    function A = asgn(A,B,I,J)
      S.type = '()';
      if nargin < 4, S.subs = {I};
      else           S.subs = {I,J};
      end
      A = subsasgn(A,S,B);
    end
    
    %--------------------------------------------------
    % Selecting sub-matrices
    %--------------------------------------------------
    function varargout = subsref(obj, S)
      switch S(1).type
        case '()'
          % Select a sub-matrix
          ind = builtin('subsref', reshape(1:numelx(obj), size(obj)), S(1) );
          out = splitvar(size(ind), [], obj.basis_(:,ind(:)));
          if length(S) > 1
            if nargout == 0
              subsref(out, S(2:end));
            else
              varargout{1:nargout} = subsref(out, S(2:end));
            end
          else
            varargout{1} = out;
          end
        case '.'
          if nargout == 0
            % Get a property, or call a member function
            % Matlab complains if you ask for an output, but the function
            % doesn't have one... so we have to pre-process in order to check
            % if there's an output or not

            % Get metaclass
            eval(sprintf('m = ?%s;',class(obj)));
            %  m = ?splitvar; 
            propInd   = find(strcmp({m.PropertyList.Name}, S(1).subs));
            methodInd = find(strcmp({m.MethodList.Name}, S(1).subs));
            if ~isempty(propInd)
              fprintf('\n'); disp(builtin('subsref', obj, S));
            elseif ~isempty(methodInd)
              % Check if the method has any outputs
              if length(m.MethodList(methodInd).OutputNames) >= 1
                % It does - print the result to the screen
                fprintf('\n'); disp(builtin('subsref', obj, S));
              else
                % It doesn't - just call the function
                builtin('subsref', obj, S);
              end
            else
              error('Unknown method or property %s', S(1).subs);
            end
          else
            [varargout{1:nargout}] = builtin('subsref', obj, S);
          end
        case '{}'
          error('Curly braces not implemented for splitvar')
      end
    end

    %--------------------------------------------------
    % Assigning sub-matrices
    %  A(S) = B
    %--------------------------------------------------
    function A = subsasgn(A, S, B)
      if isempty(A),   A = splitvar.empty(size(A));      end
      if isnumeric(A), A = splitvar(size(A), [], A(:)'); end
      switch S(1).type
        case '()'
          % Get the indices of the elements being assigned
          ind = reshape(1:numelx(A), size(A));
          ii  = subsref(ind, S(1));
          if isnumeric(B)
            A.basis_(2:end,ii) = 0;
            A.basis_(1,ii) = B(:)';
          elseif isa(B, 'splitvar')
            A.basis_ = A.basis; B.basis_ = B.basis; % Ensure same number of rows
            A.basis_(:,ii) = B.basis_;
          else
            error('Cannot assign an object of type %s to a splitvar', class(B))
          end
          
          if length(S) > 1
            A = subsref(A, S(2:end), B);
          end
        case '.'
          A = builtin('subsasgn', A, S, B);
        case '{}'
          error('Curly braces not implemented for splitvar')
      end
    end
    
    %--------------------------------------------------
    % 'end' operator
    %--------------------------------------------------
    function out = end(obj, k, n)
      assert(n <= 2, 'splitvar:end', 'Can only index up to two dimensions')
      out = obj.size_(k);
    end
    
    
    %--------------------------------------------------
    % Compute current value of the variable / equation
    %--------------------------------------------------

    % Get the index numbers of the given variables
    function varInd = getVarInd(x)
      % Confirm that this is a variable object
      assert(all(nonzeros(x.basis_) == 1) && ...
        all(sum(x.basis_) == 1) && ...
        all(x.basis_(1,:) == 0), ...
        'splitvar:set', 'Given expression is an equation, not a variable')
      
      [I,J] = find(x.basis_);
      assert(all(J == [1:size(x.basis_,2)]'), 'splitvar:set', 'Given expression is an equation, not a variable')
      varInd = I-1;
    end

    
    % Set the default value of the variable
    % Throws an error if obj is not a variable
    function set(obj, val)
      if numel(val) == 1, val = val*ones(size(obj)); end
      assert(all(size(val) == obj.size_))
      splitProb.setVarValue(obj.getVarInd, val);
      
      % Return the val via the 'val' function. If this is a symmetric
      % object, or a derived object, then the val could be different than set
%       val = obj.val;
    end
    
    % Compute the current value of the object
    function v = val(obj)     
      v = [1 splitProb.varValues'] * obj.basis;
      v = reshape(v, obj.size_);
    end    
  end
  
  methods
    %--------------------------------------------------
    % Set the function to minimize
    %--------------------------------------------------
    function minimize(f)
%       assert(isscalar(f) && f.vexity >= 0, 'splitvar:convexitychk', 'Can only minimize scalar convex functions')
      splitProb.setObjective(f);
    end
  end

  %% Testing methods
  methods
    
    function tf = issymmetric(obj)
      tf = false;
      
      % Returns true if the object is symmetric
      if obj.size_(1) ~= obj.size_(2), return; end
      if obj.size_(1) == 1, tf = true; return; end
      
      % Confirm that X' == X
      Tbasis = obj.transpose.basis;
      tf = full(max(vec(abs(Tbasis - obj.basis)))) < 1e-6;
    end
    
  end

  %%
  methods
    %--------------------------------------------------
    %--------------------------------------------------
    % Nonlinear functions
    %--------------------------------------------------
    %--------------------------------------------------
    
    %--------------------------------------------------
    % abs
    %
    % Elementwise absolute value
    %--------------------------------------------------
    function t = abs(x)
      t = zerovar(size(x,1),size(x,2));
      for i = 1:size(x,1)
        for j = 1:size(x,2)
          t = t.asgn(absFunc(x.sub(i,j)), i, j);
        end
      end
    end
    
    %--------------------------------------------------
    % max
    %
    % Syntax is the same as the built-in max function
    % t >= max(X)
    % t >= max(X, Y)
    % t >= max(X, [], DIM)
    %
    % X and Y must be affine, constant or convex
    %--------------------------------------------------
    function t = max(X, Y, DIM)
      if nargin < 3, DIM = 1; end
      if nargin < 2, Y = [];  end
      t = minmax(X, Y, DIM, @maxFunc);
    end

    %--------------------------------------------------
    % min
    %
    % Syntax is the same as the built-in min function
    % t <= min(X)
    % t <= min(X, Y)
    % t <= min(X, [], DIM)
    %
    % X and Y must be affine, constant or convex
    %--------------------------------------------------
    function t = min(X, Y, DIM)
      if nargin < 3, DIM = 1; end
      if nargin < 2, Y = [];  end
      t = minmax(X, Y, DIM, @minFunc);
    end
      
    %--------------------------------------------------
    % Internal function
    %
    % Evaluates both min and max
    %--------------------------------------------------
    function t = minmax(X, Y, DIM, func)
      if nargin < 3, DIM = 1; end
      if nargin < 2, Y = [];  end
      assert(ismember(DIM, [1 2]), 'DIM must be 1 or 2')

      if isnumeric(X), X=splitvar(size(X), [], X(:)'); end
      
      if ~isempty(Y)
        if isnumeric(Y), Y=splitvar(size(Y), [], Y(:)'); end
      
        % max(X,Y)
        assert(all(size(X) == size(Y)), 'X and Y must be the same size'); 
        t = zerovar(size(X,1),size(X,2));
        for i = 1:size(X,1)
          for j = 1:size(X,2)
            t = t.asgn(func([X.sub(i,j);Y.sub(i,j)]), i, j);
          end
        end
      else
        % max(X, [], DIM)
        if DIM == 1
          t = zerovar(1,size(X,2));
          for i = 1:size(X,2)
            t = t.asgn(func(X.sub(:,i)), 1, i);
          end
        else
          t = zerovar(size(X,1),1);
          for i = 1:size(X,1)
            t = t.asgn(func(X.sub(i,:)), i, 1);
          end
        end
      end
    end
    

    %--------------------------------------------------
    % Norms
    % 
    %  t = norm(x, p)
    %
    % t is an upper bound on the norm of x, and therefore is a convex
    % variable
    %
    % Matrix norms
    %  ... not implemented yet
    % Vector norms
    %  norm(x)      = norm(x, 2)  Euclidian norm
    %  norm(x, inf) = max(abs(x)) Infinity norm
    %  norm(x, 1)   = sum(abs(x)) 1-norm
    %
    %--------------------------------------------------
    function t = norm(x, p)
      if nargin < 2, p = []; end
      t = normFunc(x, p);
    end
    
    
    %--------------------------------------------------
    % Quadratic form
    %--------------------------------------------------
    
    
    %--------------------------------------------------
    %--------------------------------------------------
    % Sets
    %--------------------------------------------------
    %--------------------------------------------------

    % Constrain vector (x,t) to be in a second-order-cone
    %
    % ||x|| le t
    %
%     function inSecondOrderCone(x, t)
%       assert(isvector(x), 'splitvar:inSecondOrderCone', 'x must be a vector')
%       assert(isscalar(t), 'splitvar:inSecondOrderCone', 't must be a scalar')
%       assert(all(vec(x.vexity) == 0), 'splitvar:convexitychk', 'Convexity check failed; Only affine expressions can be constrained to lie in a second order cone')
%       assert(t.vexity == 0, 'splitvar:convexitychk', 'Convexity check failed; Only affine expressions can be constrained to lie in a second order cone')
%       splitProb.instance.addContainment([t;vec(x)], 'secondOrderCone');
%     end    
    
    % Constrain vector x to be in a norm ball
    %
    % ||x||_p le 1
    %
    % Where p = 1, 2 or 'inf'
%     function inNormBall(a, p)
%       if nargin < 2, p = 2; end
%       assert(isvector(a), 'splitvar:inNormBall', 'Input must be a vector - have not implemented matrix norms yet')
%       
%       assert(all(vec(a.vexity) == 0), 'splitvar:convexitychk', 'Convexity check failed; Only affine expressions can be constrained to lie in norm balls');
%       splitProb.instance.addContainment(a, 'normBall', struct('p', p));
%     end
    
%     % Constrain vector, or vectorized matrix to be non-negative
%     function inPosOrthant(a)
%       vex = vexity(a); vex = vex(:);
%       assert(~any(isnan(vex)) && all(vex <= 0), 'splitvar:convexitychk', 'Inequality constraints must be convex');
%       splitProb.instance.addContainment(vec(a), 'nonNegative');
%     end
% 
%     % Constrain symmetric matrix to be in the semidefinite cone
%     function inSemidefiniteCone(X)
%       assert(X.issymmetric, 'splitvar:inSemidefiniteCone', 'Can only constrain symmetric matrices to be in the semidefinite cone')      
% %       splitProb.instance.addContainment(symVec(X), 'semidefinite')
%       splitProb.instance.addContainment(symVec(X), 'semidefinite')
%     end
%     
%     % Constrain x to lie betwen lb and ub
%     %
%     %   lb <= x <= ub
%     %
%     % Any elements of lb and/or ub can be inf / -inf
%     %
%     function inBox(x, lb, ub)
%       assert(all(vec(x.vexity) == 0), 'splitvar:convexitychk', 'Box constraints must be affine');
%       if nargin < 3, ub = inf*ones(size(x)); end
%       assert(all(size(lb) == size(x)) && all(size(ub) == size(x)), 'splitvar:argchk', 'lb and ub must be the same size as x')
%       
%       splitProb.instance.addContainment(x(:), 'box', struct('lb', lb(:), 'ub', ub(:)));
%     end
      
    
    %--------------------------------------------------
    %--------------------------------------------------
    % Relational operators
    %--------------------------------------------------
    %--------------------------------------------------
    % All return the RHS argument, so that we can chain then -1 <= x <= 1
    function b = eq(a, b), splitProb.add(a,b,'eq'); end
    function b = ge(a, b) 
      if prod(size(a)) > 1 && issymmetric(a) || prod(size(b)) > 1 && issymmetric(b)
        splitProb.add(a, b, 'succeq');
      else
        splitProb.add(a, b, 'ge');
      end
    end
    function b = gt(a, b), b = ge(a, b); end
    function b = le(a, b)
      if prod(size(a)) > 1 && issymmetric(a) || prod(size(b)) > 1 && issymmetric(b)
        splitProb.add(a, b, 'preceq');
      else
        splitProb.add(a, b, 'le');
      end
    end
    function b = lt(a, b), b = le(a, b); end
  end

  %--------------------------------------------------
  % Convexity
  %--------------------------------------------------
  methods
    % Return a matrix of the same size as obj, where
    %   -1 == concave
    %    0 == affine
    %    1 == convex
    %  NaN == indetermined
    function vex = vexity(obj)
      % Compute the vexity of each element in the matrix
      basis = obj.basis;
      vex = sign(basis(2:end,:)) .* kron(ones(1,numelx(obj)),splitProb.getVexity);
      
      vConvex    = max(vex); % 1 if there's a convex component
      vNonConvex = min(vex); % 1 if there's a concave component

      vex = vConvex + vNonConvex;
      vex(vConvex & vNonConvex) = NaN; % Cannot determine convexity      
      
      vex = reshape(vex, obj.size_);
    end
    
    % Return the list of variables in the equation
    %  vars = getVars(x, ind)
    % vars is a sparse 0/1 matrix of size numVars * prod(size(x))
    % If ind is non-zero, returns the indices of all variables used
    % somewhere in the matrix
    function vars = getVars(x, ind)
     if nargin < 2, ind = false; end
     
     b = x.basis;
     vars = b(2:end,:) ~= 0; % Variables with non-zero multipliers
     
     if ind, vars = find(any(vars,2)); end
    end
  end

  
  %--------------------------------------------------
  % Helpers
  %--------------------------------------------------
  methods
    function disp(x)
      % Information displayed:
      % Type: affine / convex / concave / mixed convex, concave / variable

      fprintf('  <a href="matlab:help %s">%s</a>\n\n', class(x), class(x));
      
      % Check if it's a variable, an affine expression, a function or a
      % parameter
      if strcmp(class(x), 'parameter')
       cls='parameter';
      elseif ~strcmp(class(x), 'splitvar')
       cls='function'; 
      else
       B = x.basis;
       vexity = x.vexity;
       if all(B(1,:) == 0) && all(sum(abs(B)) == 1) && nnz(B) == size(B,2) && all(vexity(:) == 0)
        cls = 'variable';
       else
        cls = 'function';
       end
      end
      sz  = x.size;
      fprintf('  Split %s %ix%i', cls, sz(1), sz(2));
      
      vexity = x.vexity;
      isConvex     = vexity(:) > 0;
      isConcave    = vexity(:) < 0;
      isAffine     = vexity(:) == 0;
      isIndefinite = isnan(vexity(:));

      mixed = ' mixed';
      if all(isConvex) || all(isConcave) || all(isAffine) || all(isIndefinite)
       mixed = '';
      end
      fprintf('%s ', mixed);
      
      if ~strcmp(cls, 'variable')
       sep = '';
       if any(isConvex),     fprintf('%sconvex',     sep); sep = ' / '; end
       if any(isConcave),    fprintf('%sconcave',    sep); sep = ' / '; end
       if any(isAffine),     fprintf('%saffine',     sep); sep = ' / '; end
       if any(isIndefinite), fprintf('%sindefinite', sep); sep = ' / '; end
      end
      fprintf('\n\n');

      fprintf('  <a href="matlab:methods(''%s'')">Methods</a>\n\n', class(x));
    end
  end
end