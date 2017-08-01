#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

// custom header files
#include "../header_files/solver_cp1.h"
#include "../header_files/basic_operations.h"
#include "../header_files/projections.h"
#include "../header_files/compute_projection_conjugates.h"
#include "../header_files/compute_array_p.h"
#include "../header_files/compute_array_nu.h"
#include "../header_files/compute_array_x.h"
#include "../header_files/compute_array_xbar.h"
#include "../header_files/compute_primal_residual_chambolle.h"
#include "../header_files/compute_dual_residual_chambolle.h"

// I am going to pass everything as simple arguments,
// a structure probably would be a better solution, 
// but this one works as well and this file is not going 
// to change, i.e. one does not have to generate it

int solver_cp1(int nprox, int lr, int lc, int lir, int lic, int ar, int ac, int air, int aic, int mir, int mic, int ktr, int ktc, int ktir, int ktic, 
               int sd1r, int sd1c, int sd1ir, int sd1ic, int sd2r, int sd2c, int sd2ir, int sd2ic,  int kr, int kc, int kir, int kic, 
               int spr, int spc, int spir, int spic, int p_len, int nu_len, int x_len, int p_nu_len, int l_len, int b_len, int f_len,
               double arrayPprev[p_len], double arrayNuPrev[p_len], double arrayXprev[x_len], double arrayPnuPrev[p_nu_len],
               double arrayP[p_len], double arrayNu[nu_len],  double arrayX[x_len], double arrayPnu[p_nu_len], double arrayL[l_len], 
               double arrayXbar[x_len], double arrayXbarPrev[x_len], double arrayB[b_len], double arrayF[f_len], double matrixL[lr][lc],
               int matrixLI[lir][lic], double matrixA[ar][ac],  int matrixAI[air][aic], double matrixMI[mir][mic], double matrixK[kr][kc], 
               int matrixKI[kir][kic], double matrixKT[ktr][ktc], int matrixKTI[ktir][ktic], double matrixSD1[sd1r][sd1c], int matrixSD1I[sd1ir][sd1ic], 
               double matrixSD2[sd2r][sd2c], int matrixSD2I[sd2ir][sd2ic], double matrixSP[spr][spc], int matrixSPI[spir][spic],
               char *stringArrayProxTypes[nprox], char *stringArrayNormTypes[nprox], int arrayLenOfVectors[nprox], 
               int arrayNewProxBeg[nprox], int arrayNewProxEnd[nprox], double arrayConstants[nprox], double arrayWeights[nprox], 
               double tau, double sigma, double theta, double primaltol, double dualtol, int maxiter, double *runtime)                   
{   
     int i;
               
	 clock_t begin, end;

	 // primal residual
	 double rPrimal; 

	 // dual residual
	 double rDual; 
     
     //local K matrix I need this variable to scale the global variable hence avoiding side effects
     double matrixKcopy[kr][kc];

	 // start measure runtime 
	 begin = clock();

	 for(i = 0;i < maxiter; i++)
	 {
         // step 0: copy values
         copy_values(arrayPprev,    arrayP,    p_len);
         copy_values(arrayNuPrev,   arrayNu,   nu_len);
         copy_values(arrayXprev,    arrayX,    x_len);
         copy_values(arrayXbarPrev, arrayXbar, x_len);
         copy_values(arrayPnuPrev,  arrayPnu,  p_nu_len);
                            
         // step 1a: compute vector p: this point will be projected
         compute_array_p(lr,  lc,  matrixL, 
                         lir, lic, matrixLI,
                         l_len, arrayL, p_len, arrayP, 
                         x_len, arrayXbar, sigma);
         
         // step 1b: compute vector nu
         compute_array_nu(ar,    ac,     matrixA,
                          air,   aic,    matrixAI, 
                      	  b_len, arrayB, nu_len,  arrayNu,
                          x_len, arrayXbar, sigma);

         // step 2: compute projection         
         compute_projection_conjugates(nprox, stringArrayProxTypes, stringArrayNormTypes,
                                       arrayLenOfVectors, arrayNewProxBeg, arrayNewProxEnd,
                                       arrayConstants, arrayWeights, p_len, arrayP);
                     
         // intermediate step: create vector p_nu to simplify computation
         copy_n_values(arrayPnu, arrayP,  0, p_len,  0);               
         copy_n_values(arrayPnu, arrayNu, 0, nu_len, p_len);
                                            
         // step 3: compute array x        
         compute_array_x(mir,  mic,  matrixMI,
                         ktr,  ktc,  matrixKT,
                         ktir, ktic, matrixKTI,
                         p_nu_len, arrayPnu, x_len, arrayX, f_len, arrayF, tau); 
         
         // step 4: compute array xbar
         compute_array_xbar(x_len, arrayX, arrayXbar, arrayXprev, theta); 
         
         // intermediate step: scale the matrix by -1
         // I think I can directly export the scaled version of this matrix
         copy_matrix(kr, kc, matrixK, matrixKcopy);           
         matrix_scalar_product(kr, kc, matrixKcopy, -1.0); 

         // step 5: convergence check
         rPrimal = compute_primal_residual_chambolle(spr,  spc,  matrixSP,
                                                     spir, spic, matrixSPI,
                                                     x_len, arrayX, arrayXprev, tau);  

         
         rDual = compute_dual_residual_chambolle(kr,   kc,   matrixKcopy,  kir,   kic, matrixKI,
                                                 sd1r, sd1c, matrixSD1, sd1ir, sd1ic, matrixSD1I,
                                                 sd2r, sd2c, matrixSD2, sd2ir, sd2ic, matrixSD2I,
                                                 x_len, arrayX, arrayXbarPrev, arrayXprev,
                                                 p_nu_len, arrayPnu, arrayPnuPrev, sigma);
                                                         
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