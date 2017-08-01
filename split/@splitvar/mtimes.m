%--------------------------------------------------
% Matrix multiplication
%--------------------------------------------------
function nvar = mtimes(a, b)
if isnumeric(a), a=splitvar(size(a), [], a(:)'); end
if isnumeric(b), b=splitvar(size(b), [], b(:)'); end
if numelx(a) == 1
  nvar = a.*b;
elseif numelx(b) == 1
  nvar = a.*b;
else
  
  assert(size(a,2) == size(b,1), 'mtimes:sizechk', 'A and B matrices must have compatible dimensions')
  
  % Three cases:
  %  1. a constant
  %  2. b constant
  %  3. some a constant, some b constant
  if all(vec(a.isconstant))
    c = reshape(a.basis_(1,:),a.size_);
    b_basis = b.basis_; n = b.size_(2);
    basis = b_basis * kron(speye(n), c');
    nvar = splitvar(a.size_(1), b.size_(2), basis);
  elseif all(vec(b.isconstant))
    nvar = (b'*a')';
  else
    ca = a.isconstant; aVarCols = any(ca==0);
    cb = b.isconstant; bVarRows = any(cb==0,2)';
    
    if all(~(aVarCols & bVarRows))
      % The result is affine - so partition it into variable * const
    
      % Partition columns of a into [aConst aVariable]
      S.type = '()';
      S.subs = {':' ~aVarCols}; aConst = subsref(a, S);
      S.subs = {':'  aVarCols}; aVar   = subsref(a, S);
      
      % Partition rows of b into [bVariable; bConst]
      S.subs = { aVarCols ':'}; bConst = subsref(b, S);
      S.subs = {~aVarCols ':'}; bVar   = subsref(b, S);
      
      % a*b = aConst*bVariable + aVariable*bConst
      nvar = aConst * bVar + aVar * bConst;
    else
      % At least some of the elements are quadratic, so do this the hard way
      nvar = zerovar(size(a,1),size(b,2));
      for i = 1:size(a,1)
        for j = 1:size(b,2)
          q = dot(a.sub(i,:), b.sub(:,j));
          nvar = nvar.asgn(q, i, j);
        end
      end
    end
  end
end
