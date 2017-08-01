%--------------------------------------------------
% Elementwise right division
%--------------------------------------------------
function nvar = rdivide(a, b)
assert(isa(a, 'splitvar'), 'splitvar:rdivide', 'Cannot divide by an affine expression')
assert(numelx(a)==1 || prod(size(b))==1 || all(size(a)==size(b)), 'splitvar:times', 'Matrix dimensions must agree');
assert(isnumeric(b) || all(vec(b.isconstant)), 'splitvar:divisionchk', 'Can only divide by a constant');

if isnumeric(b)
  assert(all(vec(abs(b) > 0)), 'splitvar:dividebyzero', 'Division by zero');
else
  assert(all(vec(abs(b.basis_(1,:)) > 0)), 'splitvar:dividebyzero', 'Division by zero');
  b.basis_(1,:) = 1 ./ (b.basis_(1,:));
end
nvar = a .* (1./b);
