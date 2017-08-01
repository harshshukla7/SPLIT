#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "../header_files/projections.h"
#include "../header_files/basic_operations.h"
#include "../header_files/compute_projection.h"

// TODO: add comments

void compute_projection(int nprox, char *stringArrayProxTypes[nprox],
                        char *stringArrayNormTypes[nprox],
                        int arrayLenOfVectors[nprox], 
                        int arrayNewProxBeg[nprox], 
                        int arrayNewProxEnd[nprox], 
                        double arrayConstants[nprox], double arrayWeights[nprox],
                        int q_len, double arrayQ[q_len],
                        int y_len, double arrayY[y_len], double stepsize)
// function body
{ 
	 int i, j, length; 
     
	 // loop through all the proxes  
	 for(i = 0; i < nprox; i++) 
	 { 
		 length = arrayLenOfVectors[i];
         
		 double temp_vector[length]; 
		 double result_of_projection[length];

		 // select a chunk from vector q to project 
		 select_values(temp_vector, arrayQ, arrayNewProxBeg[i], arrayNewProxEnd[i] + 1); 

         // select the proper function to perform projection 
		 // KEEP IN MIND: return value of strcmp is zero - you have to negate it 
		 if (!strcmp(stringArrayProxTypes[i], "nonNegative")) 
			 proj_positive(result_of_projection, temp_vector, length, 1/stepsize*arrayWeights[i], stringArrayNormTypes[i], arrayConstants[i]); 
		 if (!strcmp(stringArrayProxTypes[i], "ellipse")) 
			 proj_quadball(result_of_projection, temp_vector, length, 1/stepsize*arrayWeights[i], stringArrayNormTypes[i], arrayConstants[i]);
		 if (!strcmp(stringArrayProxTypes[i], "lorentz"))
			 proj_socp(result_of_projection, temp_vector, length, 1/stepsize*arrayWeights[i], stringArrayNormTypes[i], arrayConstants[i]);
		 if (!strcmp(stringArrayProxTypes[i], "normBall"))
             proj_normBall(result_of_projection, temp_vector, length, 1/stepsize*arrayWeights[i], stringArrayNormTypes[i], arrayConstants[i]); 
 		 if (!strcmp(stringArrayProxTypes[i], "normProx"))
			 prox_norm(result_of_projection, temp_vector, length, 1/stepsize*arrayWeights[i], stringArrayNormTypes[i], arrayConstants[i]);

		 // copy n_values to y containing the solution 
		 copy_n_values(arrayY, result_of_projection, 0, length, arrayNewProxBeg[i]);
	 } 
} 