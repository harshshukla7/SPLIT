#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

// custom header files
#include "../header_files/basic_operations.h"
#include "../header_files/compute_array_p.h"

// TODO: add comments

void compute_array_p(int lr, int lc, double matrixL[lr][lc],
                     int lir, int lic, int matrixLI[lir][lic], 
                     int l_len, double arrayL[l_len],
                     int p_len, double arrayP[p_len], 
                     int x_len, double arrayXbar[x_len], double sigma)
// function body  
{
     // p = p + sigma*(L*x_bar + l)
   
     int i, j, index1, index2;
	 double sum; 
	 double Lx[lr];
	 double temp_array[l_len];
	 double temp_array1[p_len];

     // compute L*x_bar
	 for(i = 0;i < lr; i++) 
	 { 
		 index1 = matrixLI[i][0];
		 sum = 0.0; 
		 for(j = 1; j <= index1; j++) 
		 { 
			 index2 = matrixLI[i][j];
			 sum = sum + matrixL[i][index2] * arrayXbar[index2];
		 } 
		 Lx[i] = sum;
	 } 
    
     // L*x_bar + l 
     add_vectors(temp_array, Lx, arrayL, l_len);
     
     // sigma * (L*x_bar + l)
  	 scalar_vector_product(temp_array, sigma, l_len);
     
     // p + sigma * (L*x_bar + l) 
	 add_vectors(temp_array1, temp_array, arrayP, p_len);
     copy_values(arrayP, temp_array1, p_len);
}