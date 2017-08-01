function nvar = plus(a, b)

if isnumeric(a), a=splitvar(size(a), [], a(:)'); end
if isnumeric(b), b=splitvar(size(b), [], b(:)'); end
assert(numelx(a)==1 || numelx(b)==1 || all(size(a)==size(b)), 'splitvar:plus', 'Matrix dimensions must agree')
if     numelx(a)==1, a = a.*ones(size(b));
elseif numelx(b)==1, b = b.*ones(size(a));
end
nvar = splitvar(size(a), [], a.basis + b.basis);

