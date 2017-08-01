%--------------------------------------------------
% Summation - works the same as the builtin sum
%--------------------------------------------------
function nvar = sum(a, dim)
if nargin < 2, dim = 1; end
if isvector(a)
  nvar = splitvar(1, 1, sum(a.basis,2));
else
  assert(isscalar(dim) && (dim==1 || dim==2), 'splitvar:argchk', 'Second argument must be 1 or 2')
  if dim == 2
    nvar = sum(a', 1)';
  else
    a_basis = a.basis;
    nrows   = size(a,1);
    basis = spalloc(splitProb.numVars+1, 1, ceil(nnz(a_basis)/size(a_basis,2)));
    for i = 1:size(a_basis,2)/nrows
      ind = (i-1)*nrows+1:i*nrows;
      basis(:,i) = sum(a_basis(:,ind),2);
    end
    nvar  = splitvar(1, size(a,2), basis);
  end
end
