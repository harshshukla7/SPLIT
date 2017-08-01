%--------------------------------------------------
% Trace
%--------------------------------------------------
function nvar = trace(obj)
assert(obj.size_(1) == obj.size_(2), 'splitvar:trace', 'Matrix must be square')
nvar = sum(diag(obj));

