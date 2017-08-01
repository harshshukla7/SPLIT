classdef quadraticFunc < splitFunc
 %
 % q = quadraticFunc(x, Q, d, c)
 %
 %  0.5*x'*Q*x + d'*x + c = q
 %
 % Vexity of q is determined by the eigenvalues of Q
 %
 % q = quadraticFunc(x, y)
 %
 %  q = x'*y = [1;z]'*Bx'*By*[1;z]
 %
 % where z is the global variable and Bx, By are the bases of x and y
 %
 % Vexity of q is determined by the eigenvalues of Bx'*By
 %
 %
 % Calling formats:
 %
 %  quadraticFunc(x, Q, d, c)   -> q = x'*Q*x + d'*x + c
 %  quadraticFunc(x, Q, d)      -> q = x'*Q*x + d'*x
 %  quadraticFunc(x, Q)         -> q = x'*Q*x
 %  quadraticFunc(x, y)         -> q = x'*y
 %  quadraticFunc(x)            -> q = x'*x
 %
 
 properties
  % 0.5*x'*Q*x + d'*x + c = q
  x
  Q
  d
  c
  
  % sQ'*sQ = Q
  sQ
  
  % [1;z]'*H*[1;z] = q (z is all variables)
  H
 end
 
 methods
  %------------------------------------------------------------------
  % Constructor only defines inputs
  %------------------------------------------------------------------
  function func = quadraticFunc(varargin)
   func = func@splitFunc(varargin{:});
  end
  
  %------------------------------------------------------------------
  % Main initialization function
  %------------------------------------------------------------------
  function func = initialize(func, x, Q, d, c)
   if nargin < 5, c = 0; end
   if nargin < 4, d = zeros(length(x),1); end
   if nargin < 3, Q = x; end
   
   if isnumeric(x)
    % Passing in a quadratic basis
    H = x;
    [x,Q,d,c] = func.quadBasisToFunc(H);
   elseif isa(Q,'splitvar')
    assert(all(d==0) && c==0, 'd and c must both be zero for x''*y usage')
    
    % q = x'*y
    y = Q;
    
    assert(isvector(x) && isvector(y), 'x and y must be a vectors for x''*y')
    assert(length(x) == length(y), 'x and y must be column vectors')
    x = x(:); y = y(:);
    
    % Build a quadratic basis:
    %   [1;x]'*Basis_A*Basis_B'*[1;x]
    %   H := Basis_A*Basis_B'
    H = x.basis * y.basis';
    
    %   Symmetrize H := 0.5*(H+H')
    H = 0.5*(H+H');
    
    [x,Q,d,c] = func.quadBasisToFunc(H);
   end
   
   % q = 0.5*x'*Q*x + d'*x + c
   assert(size(Q,1) == size(Q,2), 'Matrix Q must be square')
   d = d(:);
   assert(size(d,1) == size(Q,1), 'Vector d must be of length %i', size(Q,1));
   assert(isscalar(c), 'c must be a scalar')
   
   % Symmetrize
   Q = sparse(Q);
   Q = (Q+Q')/2;
   
   % Build a quadratic basis:
   %   0.5*[1;z]'*Basis_x*Q*Basis_x'*[1;z] + d'*Basis_x'*[1;z] + c
   %   H := 0.5*Basis_x*Q*Basis_x' + [c d'*Basis_x'/2;Basis_x*d/2 0]
   B = x.basis;
   H = B*Q*B'/2;
   H(1,1) = H(1,1) + c;
   dd = B*d/2;
   H(:,1) = H(:,1) + dd;
   H(1,:) = H(1,:) + dd';

   % Compute the square root of Q
   % Q is semi-definite, so compute it's square root via LDL factorization
   % Compute square root via LDL factorization Q = sQ'*sQ
   [L,D,P,S] = ldl(sparse(Q));
   s  = 1./spdiags(S);
   dd = abs(sqrt(spdiags(D)));
   n  = length(dd);
   
   if n == 0
     sQ = sparse(Q);
   else
     sQ = spdiags(dd,0,n,n) * L' * P' * spdiags(s,0,n,n);
   end
   
   func.sQ = sQ;
   func.x  = x;
   func.Q  = Q;
   func.d  = d;
   func.c  = c;
   func.H  = H;
  end
  
  %------------------------------------------------------------------
  % Return the type of this function
  %------------------------------------------------------------------
  function type = getType(~)
   type = splitFuncTypes.quadratic;
  end
  
  %------------------------------------------------------------------
  % Compute the vexity of the function
  %------------------------------------------------------------------
% Removed because a linear function with a parameter in it is modeled as an
% indefinite quadratic
  function vexity = computeVexity(func)
    vexity = 1;
  end
%    Q = func.Q;
%    
%    % Check that Q is PSD or NSD
%    Q = (Q + Q')/2;
%    e = eig(Q);
%    assert(all(e >= -1e-6) || all(e <= 1e-6), 'Matrix Q must be positive semi-definite or negative semi-definite')
%    Q = sparse(Q);
%    if     all(e >= -1e-6), vexity =  1;   % Convex
%    elseif all(e <= 1e-6),  vexity = -1;   % Concave
%    else   vexity = 0;                     % Affine
%    end
%   end
  
  %------------------------------------------------------------------
  % Convert to optimizable constraint form
  %
  % One normFunc and one splitvar
  %  (Normball or SOC constraint)
  %------------------------------------------------------------------
  function con = flattenAffine(func, w, aff)
   assert(isempty(aff) && w == 1, 'Call weightedSum.simplify before flattening the function')
   
   %............
   % Ellipse
   %............
   if all(func.d == 0)
    % Ellipsoidal constraint
    %  0.5*x'*Q*x + c >= 0
    % convert to
    %  y'*y <= c
    %  y = L*x
    
    con = ellipseConSet(func.x,func.sQ,func.c*2);
   else
    % Convert to a second-order cone constraint via epigraph
    % We give no set-containment here, since the epigraph will do this for us
    con = emptyConSet(func);
    
    % [Future extension: If we have a non-zero d, then check if it can
    %  be written in the form (x-x0)'*Q*(x-x0) <= c]
   end
  end
  
  % Getter for the H
  function H = get.H(func)
   % Pad H until it matches the right number of variables
   n = splitProb.numVars;
   m = size(func.H,1);
   H = spalloc(n+1,n+1,nnz(func.H));
   H(1:m,1:m) = func.H;
  end
  
  %------------------------------------------------------------------
  % Create epigraph constraint
  %
  %  [1;x]'*H*[1;x] <= t
  %
  %------------------------------------------------------------------
  function con = flattenEpi(func)
   % Convert to an epigraph form
   %
   %   x'*Q*x + d'*x + c <= t
   % convert to soc
   %   || [1+c/2 d'/2 -1/2 ; 0 sqrtm(Q) 0] * [1;x;q] ||_2 
   %      <= [1-c/2 d'/2 1/2] * [1;x;q]
   
   n = size(func.sQ,1);
   M = [1+func.c/2  func.d'/2 -1/2; zeros(n,1) func.sQ zeros(n,1)];
   a = [1-func.c/2 -func.d'/2  1/2];

   con = socConSet(M*[1;func.x;func], a*[1;func.x;func]);
  end
  
  function quad = addLinear(quad, a)
   %------------------------------------------------------------------
   % Add a linear term to the quadratic
   %
   % quad = addLinear(quad, a)
   %
   % a is a basis for a 1x1 splitvar
   %
   %------------------------------------------------------------------
   
   H = quad.H;
   a = a(:);
   assert(size(H,1) == length(a), 'Given linear vector must be a basis')
   
   H(:,1) = H(:,1) + a/2;
   H(1,:) = H(1,:) + a'/2;
   
   quad = quad.initialize(H);
  end
 end
 
 methods (Access=private)
  function [x,Q,d,c] = quadBasisToFunc(~,H)
   % Selection basis
   ind = any(H~=0);
   ind(1) = true;
   n = size(H,1);
   
   % Reduce H to the non-zero rows/columns
   Q = H(ind,ind);
   
   % Extract the linear and constant parts
   d = 2*Q(2:end,1);
   c = Q(1,1);
   Q = Q(2:end,2:end);
   ind(1) = false;

   Q = Q * 2;
   
   P = speye(n); P = P(:,ind);
   x = splitvar(size(P,2),1,P);
  end
 end
end


