#include <stdio.h>

#include "stdio.h"
#include "string.h"
#include "splitTimer.h"
#include "splitLoad.h"
#include <Accelerate/Accelerate.h>
      
#include "splitData.h"

void test_err(double y[n]) {
  double err = 0.0;
  for(int i=0; i<n; i++) 
    err += (y[i]-sol[i])*(y[i]-sol[i]);
  if(err > 1e-8) printf("Error = %e\n", err);
}


int main(void) {
  // Load the data
  loadData();
  
  FILE *f = fopen("data.dat", "a+");

  double y[n];
  int N = 100;
  double t;
  char methods[40][40] = {"Exhaustive", "Dense", "Sparse", "Blas"};

  printf(" --- Testing n = %i, density = %.1f%% --- \n", n, density);
  for (int k=0; k<4; k++) {
    split_tic();
    if(k == 0) for(int i=0; i<N; i++) mv_exhaustive(y,x);
    if(k == 1) for(int i=0; i<N; i++) mv_dense(y,x);
    if(k == 2) for(int i=0; i<N; i++) mv_sparse(y,x);
    if(k == 3) for(int i=0; i<N; i++) mv_blas(y,x);
    t = split_toc() / (double)N;
    printf("%12s : %e seconds\n", methods[k], t);
    test_err(y);

    fprintf(f, "%s, %i, %f, %e\n", methods[k], n, density, t);
  }

  fclose(f);
  return 0;
}