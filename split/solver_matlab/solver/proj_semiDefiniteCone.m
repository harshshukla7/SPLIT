function X = proj_semiDefiniteCone(Z)
%
% X = proj_semiDefiniteCone(Z)
%
%  Project Z onto the semi-definite cone
%
%  X = argmin ||Z - X||_F s.t. X >= 0
%

isvec = false;
if isvector(Z)
  n = sqrt(length(Z));
  assert(floor(n) == ceil(n), 'proj_semiDefiniteCone:sizeChk', 'Z must be a square matrix')
  Z = reshape(Z, n, n);
  isvec = true;
end

assert(size(Z,1)==size(Z,2), 'proj_semiDefiniteCone:sizeChk', 'Z must be a square matrix')

[V,D] = eig(full(Z));
d = diag(D);

assert(all(isreal(d)), 'proj_semiDefiniteCone:symChk', 'Z must be a real symmetric matrix')

X = V*diag(max(d,0))*V';

if isvec
  X = X(:);
end
