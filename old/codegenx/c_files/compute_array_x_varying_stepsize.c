#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

// custom header files
#include "../header_files/basic_operations.h"
#include "../header_files/compute_array_x_varying_stepsize.h"
      
void compute_array_x_varying_stepsize(int ktr,   int ktc,  double matrixKT[ktr][ktc],
                                      int ktir,  int ktic, int matrixKTI[ktir][ktic], 
                                      int p_nu_len, double arrayPnu[p_nu_len], int x_len, double arrayX[x_len], 
                                      int f_len, double arrayF[f_len], int diagQ_len, double arrayDiagQ[diagQ_len], double tau)
// function body
{
     // diagQ    = diag(dat.Q)
     // diagItau = diag(eye(size(dat.Q, 1))*(1/tau))
     // inverse  = 1./(diagQ + diagItau)
     // x = inverse .* ((1/tau)*x - dat.f' - (K'*p_nu))
             
     int i, j, index1, index2;
	 double sum; 
	 double KTpn[ktr];
	 double temp_array[x_len];
	 double temp_array1[x_len];
     double temp_array2[x_len];
     double temp_array3[x_len];
     double temp_array4[x_len];

     // compute K'*p_nu
	 for(i = 0;i < ktr; i++) 
	 { 
		 index1 = matrixKTI[i][0];
		 sum = 0.0; 
		 for(j = 1; j <= index1; j++) 
		 { 
			 index2 = matrixKTI[i][j];
			 sum = sum + matrixKT[i][index2] * arrayPnu[index2];
		 } 
		 KTpn[i] = sum;
	 } 
                      
     // (1/tau)*x
     copy_values(temp_array, arrayX, x_len);
     scalar_vector_product(temp_array, 1/tau, x_len);
     
     // (1/tau)*x - f'
     subtract_vectors(temp_array1, temp_array, arrayF, x_len);
     
     // ((1/tau)*x - f' - (K'*p_nu))
     subtract_vectors(temp_array2, temp_array1, KTpn, x_len);
          
     // diagQ + diagItau 
     scalar_vector_addition(temp_array3, arrayDiagQ, 1/tau, diagQ_len);   
     
     // 1./(diagQ + diagItau)
     invert_elements(temp_array3, x_len);
     
     // inverse .* ((1/tau)*x - (K'*p_nu))        
     element_wise_vector_multiplication(temp_array4, temp_array3, temp_array2, x_len);
     copy_values(arrayX, temp_array4, x_len);
}