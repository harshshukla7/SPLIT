%--------------------------------------------------
% Right division
% Only implemented division by a constant scalar
%--------------------------------------------------
function nvar = mrdivide(a, b)
assert(isscalar(b), 'splitvar:mrdivide', 'Have only implemented division by a numerical scalar so far')
nvar = a ./ b;
