#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

// custom header files
#include "../header_files/basic_operations.h"
#include "../header_files/compute_array_nu.h"

void compute_array_nu(int ar, int ac, double matrixA[ar][ac],
                      int air, int aic, int matrixAI[air][aic], 
                      int b_len, double arrayB[b_len],
                      int nu_len, double arrayNu[nu_len], 
                      int x_len, double arrayXbar[x_len], double sigma)
// function body                      
{
     // nu = nu + sigma*(dat.A*x_bar - dat.b)
    
     int i, j, index1, index2;
	 double sum; 
	 double Ax[ar];
	 double temp_array[b_len];
	 double temp_array1[b_len];

     // compute A*x_bar
	 for(i = 0;i < ar; i++) 
	 { 
		 index1 = matrixAI[i][0];
		 sum = 0.0; 
		 for(j = 1; j <= index1; j++) 
		 { 
			 index2 = matrixAI[i][j];
			 sum = sum + matrixA[i][index2] * arrayXbar[index2];
		 } 
		 Ax[i] = sum;
	 } 
     
     // A*x_bar - b
     subtract_vectors(temp_array, Ax, arrayB, b_len);   
     
     // sigma * (A*x_bar - b)
  	 scalar_vector_product(temp_array, sigma, b_len);
     
     // nu + sigma * (A*x_bar - b) 
	 add_vectors(temp_array1, arrayNu, temp_array, b_len);
	 copy_values(arrayNu, temp_array1, nu_len);
}