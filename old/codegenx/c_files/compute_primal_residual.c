#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

// custom header files
#include "../header_files/basic_operations.h"
#include "../header_files/compute_primal_residual.h"

/*******************************************************
*                compute primal residual               *
********************************************************/

double compute_primal_residual(int slr, int slc, double matrixSL[slr][slc],
                               int slir, int slic, int matrixSLI[slir][slic], 
                               int y_len, double arrayY[y_len], double arrayPrevY[y_len],
                               int l_len, double arrayL[l_len])
{ 
	 int i, j, index1, index2;
	 double sum; 

	 // variable to store primal residual 
	 double rPrimal;

	 // array to store primal residual before computing norm 
	 double s[slr];

	 // auxiliary vector to store intermediate results
	 double temp_array[y_len];
	 double temp_array1[slr];

	 // compute y - prev_y 
	 subtract_vectors(temp_array, arrayY, arrayPrevY, y_len);

	 for(i = 0;i < slr; i++) 
	 { 
		 // store the number of nonzero elements 
		 index1 = matrixSLI[i][0];
		 sum = 0.0; 
		 for(j = 1; j <= index1; j++) 
		 { 
			 index2 = matrixSLI[i][j];
			 sum = sum + matrixSL[i][index2] * temp_array[index2];
		 } 
		 temp_array1[i] = sum;
	 } 

	 // copy primal residual into the array s 
	 copy_values(s, temp_array1, slr);

	 // compute norm 
	 compute_norm(&rPrimal, s, slr);

	 return rPrimal; 
} 