#include <stdio.h>

#include "test.h"

REAL compute_err(REAL *y, REAL *sol, int nn) {
  REAL err = 0.0;
  for(int i=0; i<nn; i++) 
    err += (y[i]-sol[i])*(y[i]-sol[i]);
  return err;
}

int main(void) {
  printf("--- Test matrix-vector multiplication methods ---\n");
  printf("This file tests for correctness of the methods, it does not compare speeds\n");

  // Load the data (which variables are loaded are in test.h)
  loadData();

  REAL y[n];

  printf("\n\n===> Testing matrix-vector products <===\n");

  func_blas(y, x);
  printf("%20s error = %.4e\n", "BLAS", compute_err(y, sol, n));

  func_forloop(y, x);
  printf("%20s error = %.4e\n", "for-loops", compute_err(y, sol, n));

  func_sparse(y, x);
  printf("%20s error = %.4e\n", "Sparse", compute_err(y, sol, n));

  func_exhaustive(y, x);
  printf("%20s error = %.4e\n", "Exhaustive", compute_err(y, sol, n));



  printf("\n\n===> Testing linear solvers <===\n");

  REAL work[x_lin_len];  // Work variable for LDL solve
  REAL y_lin[x_lin_len]; // Solution

  // We copy the input x_lin to a work variable first, since work is modified on exit
  memcpy(work, x_lin, sizeof(REAL)*x_lin_len);
  func_solve_ldl(y_lin, work);
  printf("%20s error = %.4e\n", "LDL", compute_err(y_lin, y_sol, x_lin_len));

  func_solve_inv(y_lin, x_lin);
  printf("%20s error = %.4e\n", "Invert", compute_err(y_lin, y_sol, x_lin_len));

  return 0;
}
