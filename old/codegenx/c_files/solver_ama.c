#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

// custom header files
#include "../header_files/solver_ama.h"
#include "../header_files/projections.h"
#include "../header_files/basic_operations.h"
#include "../header_files/solve_KKT_ama.h"
#include "../header_files/compute_vector_to_project.h"
#include "../header_files/compute_projection.h"
#include "../header_files/compute_over_relaxation.h"
#include "../header_files/compute_primal_residual.h"
#include "../header_files/compute_dual_residual.h"

// I am going to pass everything as simple arguments,
// a structure probably would be a better solution, 
// but this one works as well and this file is not going 
// to change, i.e. one does not have to generate it

int solver_ama(int mr, int mc, int lr, int lc, int mir, int mic, int lir, int lic, int slr, int slc, int slir, int slic, int alr, int alc, 
               int alir, int alic, int nprox, int lam_len, int t_len, int y_len, int cv_len, int x_len, int q_len, int l_len, int lx_len, 
               double matrixL[lr][lc],    int matrixLI[lr][lc], 
               double matrixM[mr][mc],    int matrixMI[mir][mic], 
               double matrixSL[slr][slc], int matrixSLI[slir][slic],  
               double matrixAlphaL[alr][alc], int matrixAlphaLI[alir][alic], 
               double arrayLam[lam_len],  double arrayT[t_len], double arrayY[y_len], double arrayPrevY[y_len], double arrayLx[lx_len],
               double arrayConstVec[cv_len], double arrayL[l_len], double arrayX[x_len], double arrayQ[q_len],
               char *stringArrayProxTypes[nprox], char *stringArrayNormTypes[nprox],
               int arrayLenOfVectors[nprox], int arrayNewProxBeg[nprox], int arrayNewProxEnd[nprox], 
               double arrayConstants[nprox], double arrayWeights[nprox],
               int maxiter, int num_opt_var, double stepsize, double dualtol, double primaltol, double alpha, double *runtime)
               
{   
     int i;
               
	 clock_t begin, end;

	 // primal residual
	 double rPrimal; 

	 // dual residual
	 double rDual; 

	 // start measure runtime 
	 begin = clock();
     
   	 for(i = 0;i < maxiter; i++)
	 {
         
         //solve linear system 
         solve_KKT_ama(mr,  mc,  matrixM, 
                       mir, mic, matrixMI,
                       t_len,    arrayT,
                       lam_len,  arrayLam,
                       cv_len,   arrayConstVec); 
                                           
         // extract NUM_OF_OPT_VAR from vector t (obtained by solving KKT) into vector x
		 select_values(arrayX, arrayT, 0, num_opt_var);
            
		 // copy y into prev_y 
		 copy_values(arrayPrevY, arrayY, y_len);   
                 
         // compute over relaxation
         compute_over_relaxation(alr,  alc,  matrixAlphaL,
                                 alir, alic, matrixAlphaLI,
                                 lx_len, arrayLx, l_len, arrayL,
                                 y_len,  arrayY,  x_len, arrayX, alpha);
         
                  
         // prepare vector we are going to project
         compute_vector_to_project(lam_len,   arrayLam,
                                   q_len,     arrayQ,
                                   l_len,     arrayL, 
                                   lx_len,    arrayLx,
                                   x_len,     arrayX, stepsize); 
    
         // compute projection
		 compute_projection(nprox,        stringArrayProxTypes, stringArrayNormTypes, arrayLenOfVectors,
                            arrayNewProxBeg, arrayNewProxEnd,         arrayConstants,        arrayWeights,
                            q_len,        arrayQ,        y_len,      arrayY,    stepsize);
         
		 // compute y - q: We need it to update lambda
		 subtract_vectors(arrayLam, arrayY, arrayQ, lam_len);
         
         // compute stepsize * lambda
		 scalar_vector_product(arrayLam, stepsize, lam_len);

		 // compute primal residual      
         /* rPrimal = compute_primal_residual(slr,    slc,    matrixSL,
                                           slir,   slic,   matrixSLI, 
                                           y_len,  arrayY, arrayPrevY,
                                           l_len,  arrayL);*/
         
         // compute dual residual
		 rDual   = compute_dual_residual(lr,     lc,      matrixL,
                                         lir,    lic,     matrixLI,
                                         y_len,  arrayY,  l_len,  arrayL,  
                                         x_len,   arrayX);

		 // if solution is within the prescribed tolerance: break 
		 // if( rDual < dualtol && rPrimal < primaltol)
         if( rDual < dualtol)
			 break;
	 } 
     
	 // calculate runtime 
	 end = clock(); 
	 *runtime = (double)(end - begin) / CLOCKS_PER_SEC;
     *runtime = 1000.0 * (*runtime);
     
     return ++i;
}