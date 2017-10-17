function gen_mv_sparse(A, Astr, dat)
% Multiply static matrix stored in compressed form and vector 
%  y = A*x

% Write data in compressed format
writeSparseMatrix(A,Astr,dat);
pl('void mv_sparse(double y[%i], const double x[%i]) {', size(A,1), size(A,2))

pl('int i=0;')
pl('for(int r=0; r<%i; r++) {', size(A,1))
pl('  y[r] = 0.0;')
pl('  for(int ic=0; ic<%s_nz_per_row[r]; ic++) {', Astr)
pl('    y[r] += %s_dat[i]*x[%s_ind[i]];', Astr, Astr)
pl('    i++;')
pl('  }')
pl('}')

pl('}')
