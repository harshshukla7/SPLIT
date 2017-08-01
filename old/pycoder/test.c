#include <stdio.h>
#include <string.h>
#include "admm.h"

int main()
{
  double par[] = {-0.1649, 0.6277};

  Sol sol;
  solve(&sol, par);

  printf("Residuals = (%g, %g)\n", sol.rPrimal, sol.rDual);

  printf("Done!\n");

  return 0;
}