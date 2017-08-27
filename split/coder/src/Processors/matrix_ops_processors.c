// Simple matrix operation library

// NOTE : Add return vectors are assumed to be allocated by caller
// NOTE : No memory violation checks are made


/*************************************************************
 ********** Things to check 
 *
 *1. vector norm 1 is wrong I think.
 *
 ************************************************************/
#include <math.h>
#include <string.h>
#include <stdio.h>
#include "user_matrix_ops.h"

static int i;

// void print_vector(const char *name, const real *x, const int len)
// {
//   printf("%s = [", name);
//   char sep = ' ';
//   forall(len) {printf("%c%g", sep, x[i]); sep=',';}
//   printf("]\n");
// }

void copy_vector( real *y, const real *x, int len){

 loop_copy_vector_matops:   for(i=0; i<len; i++){
#pragma HLS PIPELINE
        
      y[i] = x[i];  
    }
}
// returns sum(x)
real sum_vector(const real *x, const int len)
{
	real z = 0.0;
	for(int i=0; i<len; i++) z += x[i];
   return z;
}

// returns sum(x.*x)
real sum_squared(const real *x, const int len)
{
  real z = 0.0;
  forall(len) z += x[i]*x[i];
  return z;
}

// returns sum(|x|)
real sum_absolute(const real *x, const int len)
{
    real z = 0.0;
    forall(len) z += fabs(x[i]);
    return z;
}


// y = alpha*x
// Safe to scale in place (i.e., y=x)
void scale_vector(real *y, const real alpha, const real *x, const int len)
{
	for(int i=0; i<len; i++) y[i] = alpha*x[i];
}

// returns <x,y>
real dot_vector(const real *x, const real *y, const int len)
{
	real z = 0.0;
	for(int i=0; i<len; i++) z += x[i]*y[i];
   return z;
}

// add constant to each element of vector

void vector_constant_addition(real *y, const real *x, real addvalue, const int len){
	for (int i=0; i<len; i++) y[i] = x[i] + addvalue;
}
// divide each  element of vector :


void divide_vector_elements(real *y, real *x, const int len )
{

	for(int i=0; i<len; i++) y[i] = 1/x[i];
}


void scale_add_vector(real *z, real *y, real *x, real multiplayer, const int len )
{

	for(int i=0; i<len; i++) z[i] = y[i] + multiplayer*x[i];
}

// returns y = M*x
// Assumes that M is stored in row-major format
void matvec_product(const real *M, const real *x, real *y, const int nrows, const int ncols)
{
  for(int r=0; r<nrows; r++) {
    y[r] = 0.0;
    for(int c=0; c<ncols; c++) {
      y[r] += M[r*ncols + c]*x[c];
    }
  }
}


// Proximal operators

// x = argmin rho*||x||_p + 1/2*||z - x||_2^2 
void prox_norm_one(real *xprox, const real *x, const real t, const int len)
{
   
  forall(len) {
    if      (x[i] >  t) xprox[i] = x[i] - t; 
    else if (x[i] < -t) xprox[i] = x[i] + t;
    else                xprox[i] = 0.0;
  }
} 

/*void prox_norm_one(real *xprox, const real *x, const real t, const int len)
{
   
  forall(len) {
          xprox[i] = (x[i] >  t) ?  x[i] - t : 0.0;
          xprox[i] = (x[i] < -t) ?  x[i] + t : 0.0;
                    
  }
} */


 

void prox_norm_two(real *xprox, const real *x, const real t, const int len)
{ 
  real tinv, alpha, nx;
  tinv = 1.0/t;
  forall(len) xprox[i] = x[i]*tinv;
  nx = norm_two(xprox, len);
  
  alpha = 1.0 - 1.0/nx;
  
  if (nx <= 1) forall(len) xprox[i] = 0.0;
  else         forall(len) xprox[i] = x[i]*alpha;
}


/*
void prox_norm_inf(real *xprox, const real *x, const real t, const int len)
{
    real rho;
    rho = 1/t;
    forall(len) xprox[i] = x[i]*rho;
    proj_normball_one(real *xprox, real *xprox, 1, const int len);
    forall(len) xprox[i] = x[i] - t*xprox[i];
}
*/

/*void prox_norm_inf(real *xprox, const real *x, const real t, const int len)
{
    
    
}*/
/*
void prox_norm_inf(real *x, real *z, real rho)
{ t = 1/rho;
  x = z - t*proj_normBall_one(z/t, -- (used to be pdual));
}

void proj_normball_one(real *x, const real *z, const real c) 
{
  z = z/c;
 nz = norm(z, 1);
 if nz <= 1
  lam = 0;
else
  % lam solves the equation sum max(|x|-lam, 0) == 1
v = sort(abs(z));
for i = 1:length(v)
  if sum(max(v-v(i),0)) < 1, break; end
end
% We know that lam lies between v(i-1) and v(i), and that i > 1
% Solve for lam
lam = (sum(v(i:end))-1)/(length(v)-i+1);
end
x = (z+lam < 0) .* (z+lam)  + (z-lam > 0) .* (z-lam);
}
*/



/* ***************************************************************
 * Projection on norm 1 ball

// Based on :"Projection onto l1 norm ball with application to identification of sparse autoregressive model". 
   While comparing with the algorithm written note that in c all the array are shifted by 1 i.e. array starts from 
   0 not 1. */

void proj_normball_one(real *xproj, real *x, const real c, const int len)
{
    real tmp_norm, tmp_lambda, tmp_sum, tmp_glambda;
    int i;
    //int len1=len+1;
    
    
    scale_vector(xproj, 1/c, x, len);
    //printf("xproj[0] is %f and xproj[1] is %f \n", xproj[0], xproj[1]);
    tmp_norm = norm_one(xproj,len);
    
    //printf("tmp_norm is %f \n", tmp_norm);
    
    if(tmp_norm<1)
    {
        tmp_lambda=0;
       // printf("it sets lambda equal to zero. \n");
    }
    else
    {
        
        
        /* Next we find |a_k| */
        forall(len) xproj[i] = (xproj[i] < 0.0) ? -xproj[i] : xproj[i];
        
        /*We sort in ascending order */
        quick_sort (xproj, len);
        //printf("After sorting xproj[0] is %f and xproj[1] is %f \n", xproj[0], xproj[1]);
        /* Compute the table from bottom of Algorithm I */
        
        tmp_sum = -1;
        for(i = len-1; i > 0; i--){
         
            tmp_sum = tmp_sum+xproj[i];
            
           // printf("The first i is %d \n",i);
            
           // printf("temporary sum value is %f \n", tmp_sum);
            tmp_glambda = (i-len)*xproj[i-1]+tmp_sum;
            
            //printf("tmp_glambda is %f \n", tmp_glambda);
            if (tmp_glambda > 0){
                
                break;                
            }
            
        }
        
        //printf("The loop is broken at the iteration %d \n",i);
        /* if tmp_lambda >0 then oir k=i else our k=0 i.e. the sign will change for k=0 */
        if(tmp_glambda > 0)
        
        {
          //  printf("it goes with tmp_lambda is greater than zero \n");
            tmp_lambda = tmp_sum/(len-i);
        
        }
        
        else
        {   
            
          //  printf("it goes i equal to zero iteration \n");
            tmp_lambda = (tmp_sum+xproj[0])/len;
            
        }
        
        //printf("The value of the lambda is %f \n", tmp_lambda);
        /* Compute the projection using equation 8. Here I am using macro inside macro */
        forall(len) xproj[i] = (x[i]+tmp_lambda <= 0.0) ? x[i]+tmp_lambda : ((x[i]-tmp_lambda >= 0)? x[i]-tmp_lambda: 0);
    }
}



//////////////////////////////////////////
void proj_normball_two(real *xproj, const real *x, const real c, const int len)
{
  real nz = norm_two(x, len);
  //printf("c in normball two is %f \n", c);
  if (nz <= c)  copy_vector(xproj,x,len);
  else          scale_vector(xproj, c/nz, x, len);
}

void proj_normball_inf(real *xproj, const real *x, const real c, const int len)
{
  for(int i=0; i<len; i++) {
    if      (x[i] >  c) xproj[i] =  c;
    else if (x[i] < -c) xproj[i] = -c;
    else                xproj[i] = x[i];
  }
}

void proj_box(real *xproj, const real *x, const real *lb, const real *ub, const int len)
{
  for(int i=0; i<len; i++) {
    if      (x[i] > ub[i]) xproj[i] = ub[i];
    else if (x[i] < lb[i]) xproj[i] = lb[i];
    else                   xproj[i] = x[i];
  }
}

// Project the point x = [t;y] onto the second order cone
void proj_secondOrderCone(real *xproj, const real *x, const int len)
{
     // Solution from Boyd's lecture notes
  // http://www.seas.ucla.edu/~vandenbe/236C/lectures/proxop.pdf

  real t = x[0];
  //printf("t in second order cone  is %f \n", t);
  real ny = norm_two(x+1, len-1);

  if (ny <= t)       copy_vector(xproj, x, len);
  else if (ny <= -t) forall(len) xproj[i] = 0.0;
  else {
    copy_vector(xproj, x, len);
    xproj[0] = ny;
    scale_vector(xproj, (t+ny)/(2*ny), xproj, len);
  }
    
}

// Project on second-order cone conjugate
void proj_secondOrderCone_conj(real *xproj, const real *x, const int len)
{
 

  real t = -x[0];
  real ny = norm_two(x+1, len-1);

  if (ny <= t)       {copy_vector(xproj+1, x+1, len-1); xproj[0] = -x[0];}
  else if (ny <= -t) forall(len) xproj[i] = 0.0;
  else {
    copy_vector(xproj, x, len);
    xproj[0] = ny;
    scale_vector(xproj, (t+ny)/(2*ny), xproj, len);
  }
  
  xproj[0] *= -1;
}



void proj_negative(real *xproj, const real *x, const int len)
{
  forall(len) xproj[i] = (x[i] > 0.0) ? 0.0 : x[i];
}

void proj_positive(real *xproj, const real *x, const int len)
{
  forall(len) xproj[i] = (x[i] < 0.0) ? 0.0 : x[i];
}

/***********************************************************
Vector norms
 ***********************************************************/

/****** I think norm was calculated wrongly. I modified.
 *********************************************/
real norm_one(const real *x, const int len)
{
  real n = 0.0;
  for(int i=0; i<len; i++) n = x[i] > 0 ? n+x[i] : n-x[i];
    return n;
}

real norm_two(const real *x, const int len)
{ real nsqrt;
  real n = 0.0;
  for(int i=0; i<len; i++) n += x[i]*x[i];
  nsqrt = sqrt(n);
    return nsqrt;
}

real norm_inf(const real *x, const int len)
{
  real n = 0;
  for(int i=0; i<len; i++) {
    if ( x[i] > n) n = x[i];
    if (-x[i] > n) n = -x[i];
  }
  return n;
}


/**********************************************************
 *Sorting Algorithm 
 *based on: http://rosettacode.org/wiki/Sorting_algorithms/Quicksort#C
 *********************************************************/

void quick_sort (real *a, int n) {
    int i, j;
    real p, t;
    if (n < 2)
        return;
    p = a[n / 2];
    for (i = 0, j = n - 1;; i++, j--) {
        while (a[i] < p)
            i++;
        while (p < a[j])
            j--;
        if (i >= j)
            break;
        t = a[i];
        a[i] = a[j];
        a[j] = t;
    }
    quick_sort(a, i);
    quick_sort(a + i, n - i);
}


/*void prox_norm_two(real *xprox, const real *x, const real t, const int len)
{  
  real nx = norm_two(x, len);
  real alpha = 1.0 - t/nx;
  if (nx <= t) forall(len) xprox[i] = 0.0;
  else         forall(len) xprox[i] = x[i]*alpha;
} */
