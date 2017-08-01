#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

// custom header files
#include "../header_files/basic_operations.h"
#include "../header_files/compute_dual_residual_chambolle.h"

double compute_dual_residual_chambolle(int kr,    int kc,    double matrixK[kr][kc],
                                       int kir,   int kic,   int matrixKI[kir][kic],
                                       int sd1r,  int sd1c,  double matrixSD1[sd1r][sd1c],
                                       int sd1ir, int sd1ic, int matrixSD1I[sd1ir][sd1ic],
                                       int sd2r,  int sd2c,  double matrixSD2[sd2r][sd2c],
                                       int sd2ir, int sd2ic, int matrixSD2I[sd2ir][sd2ic], 
                                       int x_len, double arrayX[x_len], double arrayXbarPrev[x_len], double arrayXprev[x_len],
                                       int p_nu_len, double arrayPnu[p_nu_len], double arrayPnuPrev[p_nu_len], double sigma)
// function body                                      
{
     // r  = -K*SD1*(x_bar_prev - x_prev) + (1/sigma)*SD2*(p_nu - p_nu_prev)
    
     int i, j, index1, index2;
	 double sum; 

	 // variable to store dual residual
	 double rDual;

	 // array to store dual residual before computing norm 
	 double r[sd2r];

	 // auxiliary vector to store intermediate results
	 double temp_array[sd1r];
     double temp_array1[sd1r];
     double temp_array2[sd2r];
     double temp_array3[sd2r];
     
     // auxiliary arrays to store the solution of the subexpression
     double temp_result_expr1[kr];
     double temp_result_expr2[kr];
     
     /*********************************************************************
      *       compute expression: -K * SD1 * (x_bar_prev - x_prev)        *
      *********************************************************************/
                
     // x_bar_prev - x_prev
     subtract_vectors(temp_array, arrayXbarPrev, arrayXprev, x_len);
     
     // SD1 * (x_bar - x_prev)
	 for(i = 0;i < sd1r; i++) 
	 { 
		 // store the number of nonzero elements 
		 index1 = matrixSD1I[i][0];
		 sum = 0.0; 
		 for(j = 1; j <= index1; j++) 
		 { 
			 index2 = matrixSD1I[i][j];
			 sum = sum + matrixSD1[i][index2] * temp_array[index2];
		 } 
		 temp_array1[i] = sum;
	 } 
     
     // -K * SD1 * (x_bar - x_prev)
     matrix_vector_product(kr, kc, matrixK, temp_array1, temp_result_expr1);

    
     /********************************************************************
     *      compute expression: (1/sigma) * SD2 * (p_nu - p_nu_prev)     *
     *********************************************************************/
      
	 // p_nu - p_nu_prev
	 subtract_vectors(temp_array2, arrayPnu, arrayPnuPrev, p_nu_len);
     
     // SD2 * (p_nu - p_nu_prev)
	 for(i = 0;i < sd2r; i++) 
	 { 
		 // store the number of nonzero elements 
		 index1 = matrixSD2I[i][0];
		 sum = 0.0; 
		 for(j = 1; j <= index1; j++) 
		 { 
			 index2 = matrixSD2I[i][j];
			 sum = sum + matrixSD2[i][index2] * temp_array2[index2];
		 } 
		 temp_array3[i] = sum;
	 } 
     
     // (1/sigma) * SD2 * (p_nu - p_nu_prev)
     scalar_vector_product(temp_array3, (1/sigma), sd2r);
     copy_values(temp_result_expr2, temp_array3, sd2r);                
 
     // sub_expr1 + sub_expr2, i.e. scaledK * SD1 * (x_bar_prev - x_prev) + (1/sigma) * SD2 * (p_nu - p_nu_prev)
	 add_vectors(r, temp_result_expr1, temp_result_expr2, sd2r);

	 // compute norm 
	 compute_norm(&rDual, r, sd2r);
    
	 return rDual;   
}