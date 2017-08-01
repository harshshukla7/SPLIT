#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

// custom header files
#include "../header_files/basic_operations.h"
#include "../header_files/compute_over_relaxation.h"

/*******************************************************
*                compute over relaxation               *
********************************************************/

void compute_over_relaxation(int alr,    int alc,  double matrixAlphaL[alr][alc],
                             int alir,   int alic, int matrixAlphaLI[alir][alic], 
                             int lx_len, double arrayLx[lx_len],
                             int l_len,  double arrayL[l_len],
                             int y_len,  double arrayY[y_len], 
                             int x_len,  double arrayX[x_len], double alpha)
// function body
{
    // Lx = alphaL*x - (1-alpha)*(l-y)
    
    int i, j, index1, index2;
	double sum; 
	double temp_array[l_len];
	double temp_array1[l_len];
    double temp_array2[l_len];

    // compute alphaL*x
	for(i = 0;i < alr; i++) 
	{ 
	    index1 = matrixAlphaLI[i][0];
		sum = 0.0; 
		for(j = 1; j <= index1; j++) 
		{ 
		    index2 = matrixAlphaLI[i][j];
			sum = sum + matrixAlphaL[i][index2] * arrayX[index2];
		} 
		temp_array[i] = sum;
	 } 
     
     // l - y
     subtract_vectors(temp_array1, arrayL, arrayY, l_len);
     
     // (1-alpha)*(l - y)
     scalar_vector_product(temp_array1, 1-alpha, l_len);
     
     // alphaL*x - (1-alpha)*(l-y) 
     subtract_vectors(temp_array2, temp_array, temp_array1, l_len);
	 copy_values(arrayLx, temp_array2, lx_len);
}              