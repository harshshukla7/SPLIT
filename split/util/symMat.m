function Z = symMat(z)
%
% Z = symMat(z)
%
% Convert the vector z into a symmetric matrix Z
% i.e., undo the symetric vectorization function svec
%

n = (-1+sqrt(1+8*length(z))) / 2;
assert(floor(n) == ceil(n), 'proj_semiDefiniteCone:argchk', 'The input vector must have length N*(N+1)/2 for some N')

% Convert the input back into a symmetric matrix
Z = spalloc(n,n,nnz(z)*2);
ind = tril(ones(n))==1;
Z(ind) = z;
Z = Z + tril(Z,-1)';
