#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "../header_files/projections.h"
#include "../header_files/basic_operations.h"
#include "../header_files/compute_vector_to_project.h"

// TODO: add comments

void compute_vector_to_project(int lam_len, double arrayLam[lam_len], 
                               int q_len,   double arrayQ[q_len],
                               int l_len,   double arrayL[l_len], 
                               int lx_len,  double arrayLx[lx_len],
                               int x_len,   double arrayX[x_len], double stepsize)
// function body                           
{
     // q =  Lx + l - 1/rho * {lam, lam_hat}
         
	 double temp_lam[lam_len];
	 double temp_array[l_len];
	 double temp_array1[l_len];
     
     add_vectors(temp_array, arrayLx, arrayL, l_len); 
     copy_values(temp_lam, arrayLam, lam_len); 
     scalar_vector_product(temp_lam, 1/stepsize, lam_len);
     subtract_vectors(temp_array1, temp_array, temp_lam, l_len);
	 copy_values(arrayQ, temp_array1, q_len);
          
}    