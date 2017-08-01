#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

// custom header files
#include "../header_files/basic_operations.h"
#include "../header_files/compute_primal_residual_chambolle.h"

double compute_primal_residual_chambolle(int spr,   int spc,  double matrixSP[spr][spc],
                                         int spir,  int spic, int matrixSPI[spir][spic],
                                         int x_len, double arrayX[x_len], double arrayPrevX[x_len], double tau)
// function body   
{
     // s  = (1/tau)*(SP*(x - x_prev));
    
     int i, j, index1, index2;
	 double sum; 
    
     double rPrimal;
     double temp_array[spr];
     double temp_array1[spr];
     double s[spr];
     
     // x - x_prev
     subtract_vectors(temp_array, arrayX, arrayPrevX, x_len);
          
     // SP * (x - x_prev)
	 for(i = 0;i < spr; i++) 
	 { 
		 // store the number of nonzero elements 
		 index1 = matrixSPI[i][0];
		 sum = 0.0; 
		 for(j = 1; j <= index1; j++) 
		 { 
			 index2 = matrixSPI[i][j];
			 sum = sum + matrixSP[i][index2] * temp_array[index2];
		 } 
		 temp_array1[i] = sum;
	 } 
          
     // (1/tau) * (SP*(x - x_prev));
     scalar_vector_product(temp_array1, (1/tau), spr);
    
	 // copy primal residual into array s 
	 copy_values(s, temp_array1, spr);

	 // compute norm 
	 compute_norm(&rPrimal, s, spr);

	 return rPrimal; 
}