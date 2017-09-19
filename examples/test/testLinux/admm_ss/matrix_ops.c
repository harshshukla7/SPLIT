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
#include "matrix_ops.h"

static int i;

void print_vector(const char *name, const double *x, const int len)
{
  printf("%s = [", name);
  char sep = ' ';
  forall(len) {printf("%c%g", sep, x[i]); sep=',';}
  printf("]\n");
}


// returns sum(x)
double sum_vector(const double *x, const int len)
{
	double z = 0.0;
	for(int i=0; i<len; i++) z += x[i];
   return z;
}

// returns sum(x.*x)
double sum_squared(const double *x, const int len)
{
  double z = 0.0;
  forall(len) z += x[i]*x[i];
  return z;
}

// y = alpha*x
// Safe to scale in place (i.e., y=x)
void scale_vector(double *y, const double alpha, const double *x, const int len)
{
	for(int i=0; i<len; i++) y[i] = alpha*x[i];
}

// returns <x,y>
double dot_vector(const double *x, const double *y, const int len)
{
	double z = 0.0;
	for(int i=0; i<len; i++) z += x[i]*y[i];
   return z;
}

// returns y = M*x
// Assumes that M is stored in row-major format
void matvec_product(const double *M, const double *x, double *y, const int nrows, const int ncols)
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
void prox_norm_one(double *xprox, const double *x, const double t, const int len)
{
  forall(len) {
    if      (x[i] >  t) xprox[i] = x[i] - t;
    else if (x[i] < -t) xprox[i] = x[i] + t;
    else                xprox[i] = 0.0;
  }
}

void prox_norm_two(double *xprox, const double *x, const double t, const int len)
{
  double nx = norm_two(x, len);
  double alpha = 1.0 - t/nx;
  if (nx <= t) forall(len) xprox[i] = 0.0;
  else         forall(len) xprox[i] = x[i]*alpha;
}

/*void prox_norm_inf(double *xprox, const double *x, const double t, const int len)
{
    
    
}*/
/*
void prox_norm_inf(double *x, double *z, double rho)
{
  x = z - t*proj_normBall_one(z/t, pDual);
}

void proj_normball_one(double *x, const double *z, const double c) 
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

// Based on :Projection onto l1 norm ball with application to identification of sparse autoregressive model 
   While comparing with the algorithm written note that in c all the array are shifted by 1 i.e. array starts from 
   0 not 1. */

void proj_normball_one(double *xproj, double *x, const double c, const int len)
{
    double tmp_norm, tmp_lambda, tmp_sum, tmp_glambda;
    int i;
    //int len1=len+1;
    
    
    scale_vector(xproj, 1/c, x, len);
    printf("xproj[0] is %f and xproj[1] is %f \n", xproj[0], xproj[1]);
    tmp_norm=norm_one(xproj,len);
    
    printf("tmp_norm is %f \n", tmp_norm);
    
    if(tmp_norm<1)
    {
        tmp_lambda=0;
        printf("it sets lambda equal to zero. \n");
    }
    else
    {
        
        
        /* Next we find |a_k| */
        forall(len) xproj[i] = (xproj[i] < 0.0) ? -xproj[i] : xproj[i];
        
        /*We sort in ascending order */
        quick_sort (xproj, len);
        printf("After sorting xproj[0] is %f and xproj[1] is %f \n", xproj[0], xproj[1]);
        /* Compute the table from bottom of Algorithm I */
        
        tmp_sum=-1;
        for(i=len-1; i>0; i--){
         
            tmp_sum=tmp_sum+xproj[i];
            
            printf("The first i is %d \n",i);
            
            printf("temporary sum value is %f \n", tmp_sum);
            tmp_glambda=(i-len)*xproj[i-1]+tmp_sum;
            
            printf("tmp_glambda is %f \n", tmp_glambda);
            if (tmp_glambda>0){
                
                break;                
            }
            
        }
        
        printf("The loop is broken at the iteration %d \n",i);
        /* if tmp_lambda >0 then oir k=i else our k=0 i.e. the sign will change for k=0 */
        if(tmp_glambda>0)
        
        {
            printf("it goes with tmp_lambda is greater than zero \n");
            tmp_lambda=tmp_sum/(len-i);
        
        }
        
        else
        {   
            
          //  printf("it goes i equal to zero iteration \n");
            tmp_lambda=(tmp_sum+xproj[0])/len;
            
        }
        
        //printf("The value of the lambda is %f \n", tmp_lambda);
        /* Compute the projection using equation 8. Here I am using macro inside macro */
        forall(len) xproj[i] = (x[i]+tmp_lambda <= 0.0) ? x[i]+tmp_lambda : ((x[i]-tmp_lambda >= 0)? x[i]-tmp_lambda: 0);
    }
}



//////////////////////////////////////////
void proj_normball_two(double *xproj, const double *x, const double c, const int len)
{
  double nz = norm_two(x, len);
  if (nz <= c)  copy_vector(xproj,x,len);
  else          scale_vector(xproj, c/nz, x, len);
}

void proj_normball_inf(double *xproj, const double *x, const double c, const int len)
{
  for(int i=0; i<len; i++) {
    if      (x[i] >  c) xproj[i] =  c;
    else if (x[i] < -c) xproj[i] = -c;
    else                xproj[i] = x[i];
  }
}

void proj_box(double *xproj, const double *x, const double *lb, const double *ub, const int len)
{
  for(int i=0; i<len; i++) {
    if      (x[i] > ub[i]) xproj[i] = ub[i];
    else if (x[i] < lb[i]) xproj[i] = lb[i];
    else                   xproj[i] = x[i];
  }
}

// Project the point x = [t;y] onto the second order cone
void proj_secondOrderCone(double *xproj, const double *x, const int len)
{
  // Solution from Boyd's lecture notes
  // http://www.seas.ucla.edu/~vandenbe/236C/lectures/proxop.pdf

  double t = x[0];
  double ny = norm_two(x+1, len-1);

  if (ny <= t)       copy_vector(xproj, x, len);
  else if (ny <= -t) forall(len) xproj[i] = 0.0;
  else {
    copy_vector(xproj, x, len);
    xproj[0] = ny;
    scale_vector(xproj, (t+ny)/(2*ny), xproj, len);
  }
}

void proj_negative(double *xproj, const double *x, const int len)
{
  forall(len) xproj[i] = (x[i] > 0.0) ? 0.0 : x[i];
}

void proj_positive(double *xproj, const double *x, const int len)
{
  forall(len) xproj[i] = (x[i] < 0.0) ? 0.0 : x[i];
}

/***********************************************************
Vector norms
 ***********************************************************/

/****** I think norm was calculated wrongly. I modified.
 *********************************************/
double norm_one(const double *x, const int len)
{
  double n = 0.0;
  for(int i=0; i<len; i++) n = x[i] > 0 ? n+x[i] : n-x[i];
    return n;
}

double norm_two(const double *x, const int len)
{
  double n = 0.0;
  for(int i=0; i<len; i++) n += x[i]*x[i];
    return sqrt(n);
}

double norm_inf(const double *x, const int len)
{
  double n = 0;
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

void quick_sort (double *a, int n) {
    int i, j;
    double p, t;
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
