function Z = mpower(X,Y)
% 
% Z = mpower(X,Y)
%
% Matrix power.
%
% Valid only for Y = 2 and square matrices X
%

assert(isscalar(Y) && Y == 2, 'argchk', 'mpower is only implemented for power of two')
assert(size(X,1) == size(X,2), 'argchk', 'Matrix X must be square')

Z = X*X;
