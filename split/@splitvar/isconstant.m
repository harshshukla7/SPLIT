function tf = isconstant(obj)
%
% True if the splitvar is a constant
%
%  tf = isconstant(x)
%
% tf is a matrix of the same size as the variable
%

tf = reshape(sum(abs(obj.basis_(2:end,:))) == 0, obj.size_);
