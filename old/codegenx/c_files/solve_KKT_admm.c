#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "../header_files/projections.h"
#include "../header_files/basic_operations.h"
#include "../header_files/solve_KKT_admm.h"

/*******************************************************
*                   solve linear system                *
********************************************************/

void solve_KKT_admm(int mr, int mc, double matrixM[mr][mc],
                    int mir, int mic, int matrixMI[mir][mic], 
                    int lam_len, double arrayLam[lam_len],
                    int t_len, double arrayT[t_len], 
                    int y_len, double arrayY[y_len],
                    int cv_len, double arrayConstVec[cv_len], double stepsize)
{ 
	 int i, j, index1, index2;
	 double sum; 
     double temp_lam[lam_len];
	 double lam_y[lam_len];
	 double temp_array[cv_len];
     
	 // compute -1/rho*lam - y: ADMM
     copy_values(temp_lam, arrayLam, lam_len);
     scalar_vector_product(temp_lam, -1.0/stepsize, lam_len);
     subtract_vectors(lam_y, temp_lam, arrayY, lam_len); 
     
	 for(i = 0;i < mr; i++) 
	 { 
		 index1 = matrixMI[i][0]; 
		 sum = 0.0; 
		 for(j = 1; j <= index1; j++) 
		 { 
			 index2 = matrixMI[i][j];
             sum = sum + matrixM[i][index2] * lam_y[index2];
		 } 
		 temp_array[i] = sum;
	 } 
          
	 // add constant term to the solution obtained by solving the KKT system 
	 add_vectors(arrayT, arrayConstVec, temp_array, cv_len); 
}