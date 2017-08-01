#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

// custom header files
#include "../header_files/basic_operations.h"
#include "../header_files/compute_dual_residual.h"

/*******************************************************
*                compute dual residual                 *
********************************************************/

double compute_dual_residual(int lr, int lc, double matrixL[lr][lc],
                             int lir, int lic, int matrixLI[lir][lic], 
                             int y_len, double arrayY[y_len],
                             int l_len, double arrayL[l_len],
                             int x_len, double arrayX[x_len])
                            
{ 
	 int i, j, index1, index2;
	 double sum; 

	 // variable to store dual residual
	 double rDual;

	 // array to store dual residual before computing norm 
	 double r[l_len];

	 // auxiliary vector to store intermediate results
	 double temp_array[l_len];
	 double temp_array1[l_len];
     
	 // compute l - y
	 subtract_vectors(temp_array, arrayL, arrayY, l_len);
     
     // compute L*x
	 for(i = 0;i < lr; i++) 
	 { 
		 // store the number of nonzero elements 
		 index1 = matrixLI[i][0];
		 sum = 0.0; 
		 for(j = 1; j <= index1; j++) 
		 { 
			 index2 = matrixLI[i][j];
			 sum = sum + matrixL[i][index2] * arrayX[index2];
		 } 
		 temp_array1[i] = sum;
	 } 
     
	 // compute the final expression L*x + l - y 
	 add_vectors(r, temp_array1, temp_array, y_len);

	 // compute norm 
	 compute_norm(&rDual, r, lr);

	 return rDual; 
}