#ifndef __splitData_h__
#define __splitData_h__
#define n = 10
double x[10];
#define x_len = 10
double sol[10];
#define sol_len = 10
void mv_exhaustive(double y[10], const double x[10]) {
y[0] = 0.0;
y[1] = 0.0;
y[2] = 0.0;
y[3] = 0.0;
y[4] = 0.0;
y[5] = 0.0;
y[6] = 0.0;
y[7] = 0.0;
y[8] = 0.0 + 0.67590512805129476792*x[9];
y[9] = 0.0;
}
int A_sparse_nz_per_row[10];
#define A_sparse_nz_per_row_len = 10
int A_sparse_ind[1];
#define A_sparse_ind_len = 1
double A_sparse_dat[1];
#define A_sparse_dat_len = 1
void mv_sparse(double y[10], const double x[10]) {
int i=0;
for(int r=0; r<10; r++) {
  y[r] = 0.0;
  for(int ic=0; ic<A_sparse_nz_per_row[r]; ic++) {
    y[r] += A_sparse_dat[i]*x[A_sparse_ind[i]];
    i++;
  }
}
}
double A[100];
#define A_len = 100
void mv_dense(double y[10], const double x[10]) {
for(int r=0; r<10; r++) {
  y[r] = 0.0;
  for(int c=0; c<10; c++) y[r] += A_dense[r*10+c]*x[c];
}
}
double A_blas[100];
#define A_blas_len = 100
void mv_blas(double y[10], const double x[10]) {
cblas_dgemv(CblasRowMajor, CblasNoTrans, 10, 10, 1.0, (double*)A, 10, x, 1, 0.0, y, 1);
}


void loadData() {
 Data *dat = splitLoad("splitData.dat");
 copyVar(x, dat, "x");
 copyVar(sol, dat, "sol");
 copyVar(A_sparse_nz_per_row, dat, "A_sparse_nz_per_row");
 copyVar(A_sparse_ind, dat, "A_sparse_ind");
 copyVar(A_sparse_dat, dat, "A_sparse_dat");
 copyVar(A, dat, "A");
 copyVar(A_blas, dat, "A_blas");
freeData(dat);
}
void freeData() {
 free(x);
 free(sol);
 free(A_sparse_nz_per_row);
 free(A_sparse_ind);
 free(A_sparse_dat);
 free(A);
 free(A_blas);
}


#endif
