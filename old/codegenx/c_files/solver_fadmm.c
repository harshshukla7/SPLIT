#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

// custom header files
#include "../header_files/solver_fadmm.h"
#include "../header_files/projections.h"
#include "../header_files/basic_operations.h"
#include "../header_files/solve_KKT_fadmm.h"
#include "../header_files/compute_vector_to_project.h"
#include "../header_files/compute_projection.h"
#include "../header_files/compute_over_relaxation.h"
#include "../header_files/compute_primal_residual.h"
#include "../header_files/compute_dual_residual.h"

int solver_fadmm(int mr, int mc, int lr, int lc, int mir, int mic, int lir, int lic, int slr, int slc, int slir, int slic, int alr, int alc, 
                 int alir, int alic, int nprox, int lam_len, int t_len, int y_len, int cv_len, int x_len, int q_len, int l_len, int lx_len, 
                 double matrixL[lr][lc],    int matrixLI[lr][lc], 
                 double matrixM[mr][mc],    int matrixMI[mir][mic], 
                 double matrixSL[slr][slc], int matrixSLI[slir][slic],  
                 double matrixAlphaL[alr][alc], int matrixAlphaLI[alir][alic],
                 double arrayLam[lam_len],  double arrayPrevLam[lam_len], double arrayLamHat[lam_len],  
                 double arrayT[t_len], double arrayY[y_len], double arrayPrevY[y_len], double arrayYHat[y_len], double arrayLx[lx_len], 
                 double arrayConstVec[cv_len], double arrayL[l_len], double arrayX[x_len], double arrayQ[q_len],
                 char *stringArrayProxTypes[nprox], char *stringArrayNormTypes[nprox],
                 int arrayLenOfVectors[nprox], int arrayNewProxBeg[nprox], int arrayNewProxEnd[nprox], 
                 double arrayConstants[nprox], double arrayWeights[nprox],
                 int maxiter, int num_opt_var, double stepsize, double dualtol, double primaltol, double alpha, double *runtime)
{
     int i;

	 clock_t begin, end;

	 // parameters to decide upon Nesterov acceleration 
	 double max1, max2, Ep; 
     
     // parameters in Nesterov acceleration
	 double beta = 1.0, beta_k = 1.0;

	 // previous values of dual and primal residual
	 double rDualPrev = 1.0, rPrimalPrev = 1.0; 

	 // auxiliary arrays
	 double temp_array[y_len];
	 double temp_array1[lam_len];

	 // primal residual
	 double rPrimal = 1.0; 

	 // dual residual
	 double rDual = 1.0; 
 
	 // start measure runtime 
	 begin = clock();

	 for(i = 0;i < maxiter; i++)
	 {
		 //solve linear system             
         solve_KKT_fadmm(mr,  mc,  matrixM, 
                         mir, mic, matrixMI,
                         t_len,    arrayT,
                         lam_len,  arrayLamHat,
                         y_len,    arrayY, arrayYHat,
                         cv_len,   arrayConstVec, stepsize); 
                                           
         // extract NUM_OF_OPT_VAR from vector t (obtained by solving KKT) into vector x
		 select_values(arrayX, arrayT, 0, num_opt_var);
            
		 // copy y into prev_y 
         rDualPrev   = rDual;
         rPrimalPrev = rPrimal; 
		 copy_values(arrayPrevY, arrayY, y_len); 
         copy_values(arrayPrevLam, arrayLam, lam_len);
                 
         // compute over relaxation
         compute_over_relaxation(alr,  alc,  matrixAlphaL,
                                 alir, alic, matrixAlphaLI,
                                 lx_len, arrayLx, l_len, arrayL,
                                 y_len,  arrayY,  x_len, arrayX, alpha);
         
                  
         // prepare vector we are going to project
         compute_vector_to_project(lam_len,   arrayLamHat,
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
         rPrimal = compute_primal_residual(slr,    slc,    matrixSL,
                                           slir,   slic,   matrixSLI, 
                                           y_len,  arrayY, arrayPrevY,
                                           l_len,  arrayL);
         
         // compute dual residual
		 rDual   = compute_dual_residual(lr,     lc,      matrixL,
                                         lir,    lic,     matrixLI,
                                         y_len,  arrayY,  l_len,  arrayL,  
                                         x_len,   arrayX);
         
		 max1 = rDualPrev >= rPrimalPrev ? rDualPrev : rPrimalPrev; 
		 max2 = rDual >= rPrimal ? rDual : rPrimal; 
		 Ep = max1 - max2; 

		 // Nestrov acceleration: 
		 if(Ep > 0)
		 { 
			 beta = beta_k;
			 beta_k = (1+sqrt(4*beta*beta+1))/2;
             
			 // compute y_hat
		     subtract_vectors(temp_array, arrayY, arrayPrevY, y_len);
             scalar_vector_product(temp_array, beta - 1, y_len);
             vector_scalar_division(temp_array, beta_k, y_len);
             add_vectors(arrayYHat, arrayY, temp_array, y_len);
            
			 // compute lam_hat
			 subtract_vectors(temp_array1, arrayLam, arrayPrevLam, lam_len);
			 scalar_vector_product(temp_array1, beta - 1, lam_len);
			 vector_scalar_division(temp_array1, beta_k, lam_len);
			 add_vectors(arrayLamHat, arrayLam, temp_array1, lam_len);
		 } 
		 else 
		 { 
			 beta_k = 1;  
			 copy_values(arrayLamHat, arrayLam, lam_len);  
			 copy_values(arrayYHat, arrayY, y_len); 
		 } 

		 // if solution is within the prescribed tolerance: break 
		 if( rDual < dualtol && rPrimal < primaltol)
			 break;
	 } 

	 // calculate runtime 
	 end = clock(); 
     *runtime = (double)(end - begin) / CLOCKS_PER_SEC;
     *runtime = 1000.0 * (*runtime);
    
     return ++i;
}