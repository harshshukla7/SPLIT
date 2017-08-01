#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "../header_files/projections.h"
#include "../header_files/basic_operations.h"
#include "../header_files/solve_KKT_fama.h"

/*******************************************************
*                   solve linear system                *
********************************************************/

void solve_KKT_fama(int mr, int mc, double matrixM[mr][mc],
                    int mir, int mic, int matrixMI[mir][mic],
                    int t_len, double arrayT[t_len],   
                    int lam_len, double arrayLamHat[lam_len], 
                    int y_len, double arrayY[y_len],
                    int cv_len, double arrayConstVec[cv_len])
              
{ 
     // declare local variables
	 int i, j, index1, index2;
	 double sum; 
     double temp_lam[lam_len];
	 double lam_y[lam_len];
	 double temp_array[cv_len];
     
     // compute matrix * lam_hat (matrix: precomputed in matlab and subsequently exported to a header file)
	 for(i = 0;i < mr; i++) 
	 { 
		 index1 = matrixMI[i][0]; 
		 sum = 0.0; 
		 for(j = 1; j <= index1; j++) 
		 { 
			 index2 = matrixMI[i][j];
             sum = sum + matrixM[i][index2] * arrayLamHat[index2];
         } 
		 temp_array[i] = sum;
	 } 

	 // add constant term to the solution obtained by solving the KKT system 
	 add_vectors(arrayT, arrayConstVec, temp_array, cv_len); 
}