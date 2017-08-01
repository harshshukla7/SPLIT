function C = power(A,B)
%
% .^ Array power
%
% Implementation of standard matlab array power, but works only for squares
%
% C = A.^B
%
% B must either be 2, or a matrix of all 2's

assert(isscalar(B) && B==2 || all(size(B) == size(A)) && all(B(:) == 2), 'argchk:argumentSize', 'B must either be 2, or a matrix of all 2''s')

C = A.*A;
