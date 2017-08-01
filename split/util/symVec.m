function v = symVec(M)
%
% v = symVec(M)
%
% Return the vectorized form of the symmetric matrix M.
% Vectorizes the lower-triangular portion of the matrix.
%
% If M is of size n x n, then v will have length n*(n+1)/2
%

ind = tril(ones(size(M,1)))==1;
v = M(ind);
