function gen_mv_dense(A, Astr, dat)
% Multiply static matrix and vector using for-loops
% 
% Matrix stored in name Astr rowwise
%  y = A*x

% Store data file
dat.add('A', vec(full(A)'));

pl('void mv_dense(double y[%i], const double x[%i]) {', size(A,1), size(A,2))

pl('for(int r=0; r<%i; r++) {', size(A,1))
pl('  y[r] = 0.0;')
pl('  for(int c=0; c<%i; c++) y[r] += %s[r*%i+c]*x[c];', size(A,2), Astr, size(A,1))
pl('}')

pl('}')
