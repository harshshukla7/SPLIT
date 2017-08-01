#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "../header_files/projections.h"
#include "../header_files/basic_operations.h"
#include "../header_files/compute_projection_conjugates.h"

// TODO: add comments

void compute_projection_conjugates(int nprox, char *stringArrayProxTypes[nprox],
                                   char *stringArrayNormTypes[nprox],
                                   int arrayLenOfVectors[nprox],
                                   int arrayNewProxBeg[nprox], 
                                   int arrayNewProxEnd[nprox], 
                                   double arrayConstants[nprox], double arrayWeights[nprox],
                                   int p_len, double arrayP[p_len])
// function body
{ 
int i, j, length; 
     
	 // loop through all the proxes  
	 for(i = 0; i < nprox; i++) 
	 { 
		 length = arrayLenOfVectors[i];
         
		 double temp_vector[length]; 
		 double result_of_projection[length];

		 // select a chunk from vector p to project 
		 select_values(temp_vector, arrayP, arrayNewProxBeg[i], arrayNewProxEnd[i] + 1); 

		 // select the proper function to perform projection 
		 // KEEP IN MIND: return value of strcmp is zero - you have to negate it 
         if (!strcmp(stringArrayProxTypes[i], "nonPositive")) 
             proj_negative(result_of_projection, temp_vector, length, arrayWeights[i], stringArrayNormTypes[i], arrayConstants[i]); 
         if (!strcmp(stringArrayProxTypes[i], "ellipseConj")) 
			 prox_norm(result_of_projection, temp_vector, length, arrayWeights[i], stringArrayNormTypes[i], arrayConstants[i]);
         if (!strcmp(stringArrayProxTypes[i], "lorentzConj"))
			 proj_socp_conj(result_of_projection, temp_vector, length, arrayWeights[i], stringArrayNormTypes[i], arrayConstants[i]);
		 if (!strcmp(stringArrayProxTypes[i], "normBall"))
             proj_normBall(result_of_projection, temp_vector, length, arrayWeights[i], stringArrayNormTypes[i], arrayConstants[i]); 
 		 if (!strcmp(stringArrayProxTypes[i], "normProx"))
			 prox_norm(result_of_projection, temp_vector, length, arrayWeights[i], stringArrayNormTypes[i], arrayConstants[i]);

		 // copy n_values to y containing the solution 
		 copy_n_values(arrayP, result_of_projection, 0, length, arrayNewProxBeg[i]);
	 } 
} 