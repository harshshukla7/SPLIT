#include <stdio.h>
#include <string.h>
#include "user_ama.h"
//#ifndef __user_probData_FPGA_h__
//#include "user_probdata_FPGA.h"
//#endif
//#include "splitTimer.h"


int main()
{
 // loadData();
  
  
	// Initialize values of all variables to zero on first call
	
	// Solve the problem
	
     Sol sol;
    int number_of_solves = 1;
    Opt opt = {0.01,0.01, 20, 1};
    
    //////////////////////////// 
  //  long double time_sol_g = 0.0;
  //  long double time_total_g = 0.0;
  // long double start = split_tic();
int i;
  for ( i=0; i<number_of_solves; i++) {
    // Initialize values of all variables to zero on first call
      
	   
    initialize();
   
    solve(&sol, par_ex, &opt);
    
    
   // time_sol_g += sol.time_sol;
   // time_total_g += sol.time_total;
  }
  
  
  //long double ElapsedNano = split_toc(start)/number_of_solves; // / (double)number_of_solves;
  //sol.time_total = time_total_g/number_of_solves;
  //sol.time_sol = time_sol_g/number_of_solves;
    ////////////////////////////
    //split_tic();
    //solve(&sol, par_ex, &opt);
    //long double ElapsedNano = split_toc();

printf("Solution :\n");
printf("  Dual residual error  = %e\n", sol.rDual);
printf("  Primal residual error = %e\n", sol.rPrimal);
printf("  Numer of iterations = %d\n", (int) sol.itr);
//printf("\n\n\n Solve time is %Lf \n\n\n",ElapsedNano);
//printf("KKT solve time per iteration is %Lf \n", sol.time_sol);
//printf("total solve time per iteration is %Lf \n", sol.time_total);
//printVec(sol.primal, "Primal solution", nPrimal);
//printVec(sol.dual,   "Dual solution", nDual);
printf("number of solves are %d \n",number_of_solves);
REAL err[nPrimal];
REAL errd[nDual];

forall(nPrimal) err[i] = sol.primal[i] - sol_x[i];
printf("Error between matlab solution and c-solution for primal: %e\n", norm_two(err, nPrimal));

//forall(nDual) errd[i] = sol.aux_prim[i] - sol_y[i];
//printf("Error between matlab solution and c-solution for y: %e\n", norm_two(errd, nDual));

forall(nDual) errd[i] = sol.dual[i] - sol_lam[i];
printf("Error between matlab solution and c-solution for lambda : %e\n", norm_two(errd, nDual));

//////////////////////////////////////////////////////////////////////
////////////// Write to tex file

//FILE *f = fopen("FAMA_final.txt", "a+");
//if (f == NULL)
//{
//    printf("Error opening file!\n");
//    exit(1);
//}

/* print integers and floats */
//fprintf(f, "%d %d %d %d %Lf %Lf %Lf %Lf %f %f \n", nPrimal, nDual, number_of_solves, sol.itr , sol.time_sol, sol.time_total, ElapsedNano, (sol.time_total-sol.time_sol), norm_two(err, nPrimal), norm_two(errd, nDual));
//fprintf(f, "nPrimal is: %d, nDual is: %d and time in nanosecond is : %Lf\n", nPrimal, nDual, ElapsedNano);

//fclose(f);

//////////////////////////////////////////////////////////////////////
////////////// Write to tex file ends

  return 0;
}
