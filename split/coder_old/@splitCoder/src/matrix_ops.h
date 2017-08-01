// Simple matrix operation library

// NOTE : Add return vectors are assumed to be allocated by caller
// NOTE : No memory violation checks are made

#ifndef __matrix_ops__
#define __matrix_ops__
#include <math.h>
#include <string.h>

void print_vector(const char *name, const double *x, const int len);

// y = x
#define copy_vector(y,x,len) memcpy(y,x,sizeof(double)*len)

// returns sum(x)
double sum_vector(const double *x, const int len);

// returns sum(x.*x)
double sum_squared(const double *x, const int len);

// y = alpha*x
void scale_vector(double *y, const double alpha, const double *x, const int len);

// returns <x,y>
double dot_vector(const double *x, const double *y, const int len);

// y = M*x
void matvec_product(const double *M, const double *x, double *y, const int nrows, const int ncols);

// Macro to extend an operation elementwise
#define forall(range_h) for(i=0; i<range_h; i++)

// Proximal functions
void prox_norm_one(double *xprox, const double *x, const double t, const int len);
void prox_norm_two(double *xprox, const double *x, const double t, const int len);
//void prox_norm_inf(double *x, double *z, double rho);
//void proj_normball_one(double *x, const double *z, const double c);

// Projections
void proj_normball_two(double *xproj, const double *x, const double c, const int len);
void proj_normball_inf(double *xproj, const double *x, const double c, const int len);
void proj_secondOrderCone(double *xproj, const double *x, const int len); // Project the point x = [t;y] onto the second order cone
void proj_negative(double *xproj, const double *x, const int len);
void proj_positive(double *xproj, const double *x, const int len);

// Not used
void proj_box(double *xproj, const double *x, const double *lb, const double *ub, const int len);

// Vector norms
double norm_one(const double *x, const int len);
double norm_two(const double *x, const int len);
double norm_inf(const double *x, const int len);

#endif
