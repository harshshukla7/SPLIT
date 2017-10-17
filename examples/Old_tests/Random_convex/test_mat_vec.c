#include <stdio.h>
#include <string.h>
#include "matrix_ops.h"

#include <Accelerate/Accelerate.h>
#include "splitTimer.h"
#include "probData.h"

int main(){

loadData();
int i;
long double start_time, end_time, sum_time;
sum_time = 0.0;

start_time = split_tic();
for (i=0;i<itr;i++){
   
custom_mult_A(vec_c,vec_b);

custom_mult_A(vec_b,vec_c);

}
end_time = split_toc(start_time);



double resi =0.0;
forall(vec_b_len) resi += fabs(vec_b[i]-sol[i]);

//printf("residual is %f and time is %Lf \n",resi, end_time);



//////////////////////////////////////////////////////////////////////
////////////// Write to tex file

FILE *f = fopen("MAT_VEC_text_exchustive.txt", "a+");
if (f == NULL)
{
    printf("Error opening file!\n");
    exit(1);
}

/* print integers and floats */
fprintf(f, "%d %d %d  %d %d %d %Lf %Lf %Lf %f \n", vec_b_len, bd, itr, zerosA, nonzerosA, totalA, end_time/(2*itr), end_time/(itr), end_time, resi);
//fprintf(f, "nPrimal is: %d, nDual is: %d and time in nanosecond is : %Lf\n", nPrimal, nDual, ElapsedNano);

fclose(f);


return 0;

}

//gcc -std=c99 -Wall -o mat test_mat_vec.c splitLoad.c probData.c