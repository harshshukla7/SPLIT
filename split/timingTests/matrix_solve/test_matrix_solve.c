#include <stdio.h>

#include "splitTimer.h"
#include "splitLoad.h"
#include <Accelerate/Accelerate.h>

#include "splitData.h"

void test_err(double y[n]) {
  double err = 0.0;
  for(int i=0; i<n; i++) 
    err += (y[i]-sol[i])*(y[i]-sol[i]);
  if(err > 1e-8) printf("Error = %e\n", err);
  printf("Err = %e\n", err);
  // printVec(y, "y", n);
  // printVec(sol,"sol", n);
}

int main() {
  FILE *f = fopen("data.dat", "a+");

  loadData();

  double x[b_len];
  double b0[b_len];

  memcpy(b0,b,sizeof(double)*b_len);

  int N = 100;
  double t;
  char methods[40][40] = {"Invert", "LDL", "Banded"};

  printf(" --- Testing n = %i, density = %.1f%%, horizon = %i --- \n", n, density, horizon);
  for (int k=0; k<3; k++) {
    split_tic();
    if(k == 0) for(int i=0; i<N; i++) solveKKT_invert(x,  b);  
    if(k == 1) for(int i=0; i<N; i++) {memcpy(b,b0,sizeof(double)*b_len); solveKKT_ldl(x, b);}
    if(k == 2) for(int i=0; i<N; i++) {solveKKT_banded(x, b);}
    t = split_toc() / (double)N;
    printf("%12s : %e seconds\n", methods[k], t);
    test_err(x);

    fprintf(f, "%s, %i, %f, %i, %e\n", methods[k], n, density, horizon, t);
  }

  fclose(f);
  return 0;
}