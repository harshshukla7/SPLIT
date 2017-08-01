%--------------------------------------------------
% Transpose - ctranspose
%--------------------------------------------------
function nvar = ctranspose(a)
ind   = reshape(1:numelx(a), size(a))';
b     = a.basis;
nvar  = splitvar(fliplr(size(a)), [], b(:,ind(:)));
