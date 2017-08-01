%--------------------------------------------------
% Diff - works the same as builtin diff
%--------------------------------------------------
function nvar = diff(obj, n, dim)
if nargin < 3, dim = 1; end
if nargin < 2, n = 1; end

basis     = sparse(1,1);
obj_basis = obj.basis_;
sz        = obj.size_;
for i = 1:size(obj_basis,1)
  x = reshape(obj_basis(i,:), sz);
  x = diff(x, n, dim);
  if i == 1, basis = sparse(size(obj_basis,1), numel(x)); end
  basis(i,:) = x(:)';
end
nvar = splitvar(size(x), [], basis);
