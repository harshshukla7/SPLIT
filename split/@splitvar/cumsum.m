%--------------------------------------------------
% Cumulative sum - works the same as builtin cumsum
%--------------------------------------------------
function nvar = cumsum(obj, dim)
if isvector(obj)
  nvar = splitvar(obj.size_, [], cumsum(obj.basis_, 2));
else
  if nargin < 2, dim = 1; end
  basis = obj.basis_; obj_basis = obj.basis_; sz = obj.size_;
  for i = 1:size(obj_basis,1)
    x = reshape(obj_basis(i,:), sz);
    x = cumsum(x, dim);
    basis(i,:) = x(:)';
  end
  nvar = splitvar(sz, [], basis);
end
