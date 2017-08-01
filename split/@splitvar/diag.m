%--------------------------------------------------
% Diagonal - works the same as builtin diag
%--------------------------------------------------
function nvar = diag(obj, K)
if isvector(obj)
  n = numelx(obj);
  basis = sparse(splitProb.numVars+1, n*n);
  ind = [1:n:n*n] + [0:n-1]; % Indices of diagonal elements
  basis(:,ind) = obj.basis;
  nvar = splitvar(n, n, basis);
else
  if nargin < 2, K = 0; end
  ind = reshape(1:numelx(obj), obj.size_);
  basis = obj.basis_(:,diag(ind,K));
  nvar = splitvar(size(basis,2), 1, basis);
end
