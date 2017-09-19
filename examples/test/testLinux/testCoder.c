#include <stdio.h>
#include <string.h>
#include "admm.h"

int main()
{
  loadData();

	// Initialize values of all variables to zero on first call
	initialize();

	// Solve the problem
	Sol sol;
	Opt opt = {1e-4, 1e-4, 1000, 10};
	solve(&sol, par_ex, &opt);

printf("Solution :\n");
printf("  Dual residual   = %e\n", sol.rDual);
printf("  Primal residual = %e\n", sol.rPrimal);
printf("  Numer of iterations = %i\n", sol.itr);
printVec(sol.primal, "Primal solution", nPrimal);
printVec(sol.dual,   "Dual solution", nDual);

REAL err[nPrimal];
int i;
forall(nPrimal) err[i] = sol.primal[i] - sol_x[i];
printf("Error between matlab solution and c-solution : %e\n", norm_two(err, nPrimal));

  return 0;
}
