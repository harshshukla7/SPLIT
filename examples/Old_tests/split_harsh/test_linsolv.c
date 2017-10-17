#include <stdio.h>
#include <string.h>
#include "matrix_ops.h"

#include <Accelerate/Accelerate.h>

#include "probData.h"
#include "ldl.h"
int main(){

static double x[3]; // Extra rows are working space for solving KKT system
loadData();


#ifdef ldl
custom_compute_prefactor();

custom_KKT_solve(x, rhs);

#endif

#ifdef chol



#ifdef suitesparse_linsolve
custom_compute_prefactor();
#endif

custom_chol_solve(x, rhs);

printf(" rhs is %f, %f and %f \n", rhs[0], rhs[1],rhs[2]);
printf(" x is %f, %f and %f \n", x[0], x[1],x[2]);

//printf(" Lp is %d, %d and %d \n", Lp_chol[0], Lp_chol[1],Lp_chol[2]);
//printf(" Li is %d, %d and %d \n", Li_chol[0], Li_chol[1], Li_chol[2]);
//printf(" Lx is %f, %f and %f \n", Lx_chol[0], Lx_chol[1], Lx_chol[2]);

#endif

return 0;

}