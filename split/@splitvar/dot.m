function y = dot(a,b)
%
% y = dot(a,b)
%
% Compute the dot product of a and b
%
% Equivalent to a'*b
%

if isnumeric(a), a=splitvar(size(a), [], a(:)'); end
if isnumeric(b), b=splitvar(size(b), [], b(:)'); end

assert(isvector(a) && isvector(b), 'dot:sizechk', 'a and b must be vectors')
a = a(:); b = b(:);
assert(length(a) == length(b), 'dot:sizechk', 'a and b must be the same length')

if all(a.isconstant) || all(b.isconstant)
  % Linear term
  y = a'*b;
else
  % Quadratic term
  y = quadraticFunc(a, b);
end
