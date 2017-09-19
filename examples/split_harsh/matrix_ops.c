// Simple matrix operation library

// NOTE : Add return vectors are assumed to be allocated by caller
// NOTE : No memory violation checks are made

#include <math.h>
#include <string.h>
#include <stdio.h>
#include "matrix_ops.h"
    
#include "probData.h"

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
double norm_one(const double *x, const int len)
{
  double n = 0.0;
  for(int i=0; i<len; i++) n = x[i] > 0 ? x[i] : -x[i];
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