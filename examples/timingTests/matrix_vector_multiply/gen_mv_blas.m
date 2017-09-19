function gen_mv_dense(A, Astr, dat)
% Multiply static matrix and vector using BLAS


% Store data file
dat.add(Astr, vec(full(A)'));

pl('void mv_blas(double y[%i], const double x[%i]) {', size(A,1), size(A,2))

pl('cblas_dgemv(CblasRowMajor, CblasNoTrans, %i, %i, 1.0, (double*)A, %i, x, 1, 0.0, y, 1);',...
  size(A,1), size(A,2), size(A,1));

pl('}')
