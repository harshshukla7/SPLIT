function nvar = times(a, b)

% Elementwise addition

if isnumeric(a), a=splitvar(size(a), [], a(:)'); end
if isnumeric(b), b=splitvar(size(b), [], b(:)'); end
assert(numelx(a)==1 || numelx(b)==1 || all(size(a)==size(b)), 'splitvar:times', 'Matrix dimensions must agree');

if all(vec(a.isconstant | b.isconstant))
  
  % Each variable must be the product of one constant, and one variable / constant
  ba = a.basis;
  bb = b.basis;
  if     numelx(a) == 1
    if a.isconstant, basis = bb * ba(1);
    else             basis = ba * bb(1,:);
    end
  elseif numelx(b) == 1
    if b.isconstant, basis = ba * bb(1);
    else             basis = bb * ba(1,:);
    end
  else
    a_basis = a.basis; b_basis = b.basis;
    basis = sparse(size(a_basis,1),size(a_basis,2));
    
    caInd = sum(abs(a_basis(2:end,:))) == 0;
    if any(caInd)
      ca    = a_basis(1,caInd);
      bb    = b_basis(:,caInd) * diag(sparse(ca));
      basis(:,caInd) = bb;
    end
    
    cbInd = ~caInd;
    if any(cbInd)
      cb    = b_basis(1,cbInd);
      ba    = a_basis(:,cbInd) * diag(sparse(cb));
      basis(:,cbInd) = ba;
    end
  end
  
  nvar = splitvar(max([size(a);size(b)]), [], basis);
  
else
  % Quadratic form a*b for a and b scalar
  nvar = zerovar(size(a,1),size(b,2));
  for i = 1:size(a,1)
    for j = 1:size(b,2)
      nvar = nvar.asgn(dot(a.sub(i,j), b.sub(i,j)), i, j);
    end
  end
end
